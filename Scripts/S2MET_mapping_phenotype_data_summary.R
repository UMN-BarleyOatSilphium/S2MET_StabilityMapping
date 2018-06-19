## S2MET Mapping
## 
## Phenotypic data summary and analysis
## This notebook will provide some phenotyping data summaries for the S2MET Mapping project. It will include:

# 1. Basic model for g + e + gxe
# 2. Heritability estimates
# 3. Correlations among environments

# Repository directory
repo_dir <- getwd()
# Project and other directories
source(file.path(repo_dir, "source.R"))

# Other packages
library(lme4)
library(ggridges)
library(modelr)
library(pbr)


### Basic Summaries
### 
### Look at the number of lines per environment

# Find the total number of possible line x environment combinations and find
# the proportion that are observed for each trait
observed <- S2_MET_BLUEs_use %>% 
  distinct(trait, environment, line_name) %>%
  mutate(observed = TRUE) %>%
  group_by(trait) %>%
  complete(environment, line_name, fill = list(observed = FALSE)) 

(prob_observed <- observed %>%
  # group_by(trait) %>%
  summarize(prop_obs = mean(observed)))

## For each line, calculate the balance - sort lowest to highest
line_observed <- observed %>% 
  group_by(trait, line_name) %>% 
  summarize(p_env = mean(observed)) %>% 
  arrange(p_env)
# How many lines were observed everywhere?
line_observed %>% filter(p_env == 1) %>% summarize(n_line = n_distinct(line_name))

## For each environment, calculate the balance
env_observed <- observed %>% 
  group_by(trait, environment) %>% 
  summarize(p_line = mean(observed)) %>% 
  arrange(p_line)


## Number of total environments and number of environments per trait
n_distinct(S2_MET_BLUEs_use$environment)

S2_MET_BLUEs_use %>%
  group_by(trait) %>%
  summarize(n_env = n_distinct(environment))

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

# Plot
(g_met_dist <- S2_MET_BLUEs_use %>%
  mutate(environment = parse_factor(environment, levels = env_order)) %>%
  ggplot(aes(x = value, y = environment, fill = environment)) +
  geom_density_ridges() +
  facet_grid(. ~ trait, scales = "free_x") +
  scale_fill_discrete(guide = FALSE) +
  ylab("Environment") +
  xlab("") +
  labs(title = "Trait Distributions in All Environments") +
  theme_bw() )

# Sort on latitute
env_order <- S2_MET_BLUEs_use %>% 
  left_join(trial_info) %>% 
  distinct(environment, latitude) %>%
  arrange(latitude, environment) %>%
  # distinct(environment, longitude) %>% 
  # arrange(longitude, environment) %>% 
  pull(environment)

(g_met_dist <- S2_MET_BLUEs_use %>%
    mutate(environment = parse_factor(environment, levels = env_order)) %>%
    ggplot(aes(x = value, y = environment, fill = environment)) +
    geom_density_ridges() +
    facet_grid(. ~ trait, scales = "free_x") +
    scale_fill_discrete(guide = FALSE) +
    ylab("Environment") +
    xlab("") +
    labs(title = "Trait Distributions in All Environments") +
    theme_bw() )


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
  do(suppressWarnings(herit(object = .$full_fit[[1]], n_e = .$n_e, n_r = .$n_r,
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
       caption = "Estimates have been corrected for bias and error bars reflect\na 95% confidence interval of 100 bootstrap replications.") +
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

# Display nicer
prop_varcomp %>% 
  select(-vcov) %>% 
  spread(grp, prop_vcov)


# Plot
g_prop_varcomp <- prop_varcomp %>%
  ggplot(aes(x = trait, y = prop_vcov, fill = grp)) +
  geom_col() +
  ylab("Proportion") +
  xlab("Trait") + 
  scale_fill_discrete(guide = guide_legend(title = "Variance\nComponent")) +
  labs(title = "Components of Phenotypic Variance") +
  theme_bw()


ggsave(filename = "var_comp_prop.jpg", plot = g_prop_varcomp, path = fig_dir,
       height = 5, width = 5, dpi = 1000)

## Output a table
prop_varcomp_table <- prop_varcomp %>% 
  mutate(vcov_prop = str_c(round(vcov, 2), " (", round(prop_vcov * 100, 2), "%)")) %>% 
  select(-vcov:-prop_vcov) %>% 
  spread(trait, vcov_prop)

write_csv(x = prop_varcomp_table, path = file.path(fig_dir, "trait_varcomp.csv"))


### Calculate the proportion of GxE that is due to environmental genetic variance
### heterogeneity versus lack of environmental correlation

# For each environment, calculate the genetic variance via reml
env_varG <- S2_MET_BLUEs_use %>% 
  group_by(trait, environment) %>%
  do(varG = {
    df <- .
    
    lmer_control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
    wts <- df$std_error^2
  
    fit <- lmer(formula = value ~ 1 + (1|line_name), data = df, control = lmer_control,
                weights = wts)  
    
    as.data.frame(VarCorr(fit))[1,"vcov"]
    
  })

# Calculate the heterogeneity
env_varG_V <- env_varG %>% 
  group_by(trait) %>% 
  unnest() %>% 
  summarize(V = var(sqrt(varG)))

# Use the estimate of varGE across all environments to calculate L
env_L <- left_join(env_varG_V, subset(prop_varcomp, grp == "GxE",c(trait, vcov))) %>% 
  mutate(L = vcov - V)

# Use the estimate of genetic variance across all environments to calculate the 
# genetic correlation
env_r <- left_join(env_L, subset(prop_varcomp, grp == "Genotype", c(trait, vcov)), by = "trait") %>% 
  mutate(r_G = vcov.y / (vcov.y + L))

## Genetic correlation across environments:
## GY: 0.203
## HD: 0.726
## PH: 0.334


## What proportion do V and L make up of varGE?
## This is from Li et al 2018
env_r %>% 
  select(trait, varGE = vcov.x, V, L) %>% 
  mutate_at(vars(V, L), funs(prop = . / varGE))






