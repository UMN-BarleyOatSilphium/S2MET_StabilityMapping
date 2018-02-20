## S2MET Mapping
## 
## Phenotypic data summary and analysis
## This notebook will provide some phenotyping data summaries for the S2MET Mapping project. It will include:

# 1. Basic model for g + e + gxe
# 2. Heritability estimates
# 3. Correlations among environments

library(tidyverse)
library(broom)
library(stringr)
library(readxl)
library(modelr)
library(pbr)
library(rrBLUP)
library(ggridges)

# Project and other directories
source("C:/Users/Jeff/Google Drive/Barley Lab/Projects/S2MET_Mapping/source.R")

tp_geno <- tp_geno_multi

load(file.path(geno_dir, "S2_genos_hmp.RData"))

# Filter environments for those in which the TP was observed
S2_MET_BLUEs_use <- S2_MET_BLUEs %>% 
  filter(line_name %in% tp_geno,
         trait != "TestWeight")


### Basic Summaries
### 
### Look at the number of lines per environment


# Find the total number of possible line x environment combinations and find
# the proportion that are observed for each trait
(prob_observed <- S2_MET_BLUEs_use %>% 
  distinct(trait, environment, line_name) %>%
  mutate(observed = TRUE) %>% 
  group_by(trait) %>%
  complete(environment, line_name, fill = list(observed = FALSE)) %>%
  # group_by(trait) %>%
  summarize(prop_obs = mean(observed)))


## Interpretation: of the possible line x environment combinations in which at least one line was observed, the above is the proportion of the lines that were observed.



## Visualization of distributions

# Sort on grain yield environmental mean
env_order <- S2_MET_BLUEs_use %>% 
  mutate_if(is.character, as.factor) %>% 
  group_by(environment, trait) %>% 
  mutate(env_mean = mean(value, na.rm = TRUE)) %>% 
  filter(trait == "GrainYield") %>% 
  complete(environment) %>%
  arrange(env_mean) %>%
  pull(environment) %>% 
  unique()

(g_met_dist <- S2_MET_BLUEs_use %>%
  filter(trait != "TestWeight") %>%
  mutate(environment = parse_factor(environment, levels = env_order)) %>%
  ggplot(aes(x = value, y = environment, fill = environment)) +
  geom_density_ridges() +
  facet_grid(. ~ trait, scales = "free_x") +
  scale_fill_discrete(guide = FALSE) +
  ylab("Environment") +
  xlab("") +
  labs(title = "Trait Distributions in All Environments") +
  theme_bw() )

# Sort

# Save it
save_file <- file.path(fig_dir, "met_trait_dist.jpg")
ggsave(filename = save_file, plot = g_met_dist, width = 7, height = 5)



### Heritability


# Use a different optimizer
library(optimx)

# Group by trait and fit the multi-environment model
stage_two_model_fits <- S2_MET_BLUEs_use %>% 
  group_by(trait) %>%
  # mutate(value = scale(value)) %>% # Scale and center the observations
  rename(env = environment) %>% # Rename environment
  do({
    # Extract the df
    df <- .
    
    # List of models
    models <- formulas(~ value,
                      base = ~ 1 + (1|env) + (1|line_name),
                      full = add_predictors(base, ~ (1 | line_name:env)))
    
    # Lmer control
    lmer_control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore",
                                optimizer ='optimx', optCtrl=list(method='nlminb'))
    # Get the weights
    wts <- df$std_error^2
    
    # Fit the models
    fits <- models %>% 
      map(~lmer(formula = ., data = df, control = lmer_control, weights = wts))
    
    # Table of lines by environments (i.e. plots)
    plot_table <- xtabs(formula = ~ line_name + env, data = df)
    
    # Find the number of environments
    n_e <- plot_table %>%
      ifelse(. > 1, 1, .) %>%
      rowSums() %>% 
      harm_mean()
    
    # Now replicates
    n_r <- plot_table %>% 
      harm_mean()
    
    # Return data_frame
    data_frame(fits = list(fits), n_e = n_e, n_r = n_r) })
    


stage_two_lrt <- stage_two_model_fits %>% 
  mutate(llik = map(fits, ~map_df(., logLik))) %>% 
  unnest(llik) %>%
  select(-fits:-n_r) %>%
  mutate(lr = -2 * (base - full),
         p_value = pchisq(q = lr, df = 1, lower.tail = FALSE) / 2)



# The interaction term is significant for all traits
# The broad-sense heritability is calculated on a entry-mean basis according to the formula

# Now calculate heritability across all environments

# Extract the full model and calculate heritability
stage_two_herit <- stage_two_model_fits %>% 
  mutate(full_fit = map(fits, "full")) %>%
  do(suppressWarnings(herit_boot(object = .$full_fit[[1]], n_e = .$n_e, n_r = .$n_r, boot.reps = 10,
                       exp = "line_name / (line_name + (line_name:env / n_e) + (Residual / n_r))")))
    

# Plot
g_herit <- stage_two_herit %>% 
  mutate(herit_correct = heritability - bias) %>%
  ggplot(aes(x = trait, y = herit_correct, ymin = ci_lower, ymax = ci_upper)) +
  geom_col(aes(fill = "fill"), width = 0.80) +
  geom_errorbar(width = 0.5) +
  scale_fill_discrete(guide = FALSE) +
  ylab("Heritability") +
  xlab("Trait") +
  labs(title = "Broad-Sense Heritability",
       caption = "Estimates have been corrected for bias and error bars reflect\na 95% confidence interval of 500 bootstrap replications.") +
  theme_bw()

save_file <- file.path(fig_dir, "heritability.jpg")
ggsave(filename = save_file, plot = g_herit, height = 5, width = 4)


# What is the proportion of each variance component to the total phenotypic variance?

prop_varcomp <- stage_two_model_fits %>% 
  mutate(full_fit = map(fits, "full")) %>% 
  do(var_comp = as.data.frame(VarCorr(.$full_fit[[1]]))) %>% 
  unnest() %>%
  select(trait:grp, vcov) %>% 
  group_by(trait) %>%
  mutate(grp = str_replace_all(grp, c("line_name:env" = "GxE", "line_name" = "Genotype", "env" = "Environment")),
         grp = factor(grp, levels = c("Genotype", "GxE", "Environment", "Residual")),
         prop_vcov = vcov / sum(vcov))

# Plot
g_prop_varcomp <- prop_varcomp %>%
  ggplot(aes(x = trait, y = prop_vcov, fill = grp)) +
  geom_col() +
  ylab("Proportion") +
  xlab("Trait") + 
  scale_fill_discrete(guide = guide_legend(title = "Variance\nComponent")) +
  labs(title = "Components of Phenotypic Variance") +
  theme_bw()


save_file <- file.path(fig_dir, "var_comp_prop.jpg")
ggsave(filename = save_file, plot = g_prop_varcomp, height = 5, width = 5)



