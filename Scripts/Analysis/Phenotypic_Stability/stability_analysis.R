## S2MET Mapping
## 
## Phenotypic stability analysis

library(tidyverse)
library(readxl)
library(lme4)
library(broom)
library(stringr)
library(modelr)
# library(FW)
library(ggridges)
library(ggforce)
library(pbr)
library(qvalue)
library(cowplot)
library(patchwork)
library(boot)
library(rrBLUP)
library(lme4qtl)

# Repository directory
repo_dir <- getwd()

# Project and other directories
source(file.path(repo_dir, "source.R"))



# Load fw data
load(file.path(result_dir, "S2MET_pheno_mean_fw_results.RData"))



### Finlay-Wilkinson Regression

## Results

# What are the ranges of genotype values and environmental effects per trait?

pheno_summary <- S2_MET_pheno_mean_fw %>% 
  distinct(trait, line_name, environment, h, g)  %>% 
  gather(term, value, g:h) %>% 
  group_by(trait, term) %>% 
  summarize_at(vars(value), funs(min, max)) %>% 
  mutate(range = max - min)

## Environment metadata
trial_info %>% 
  filter(environment %in% S2_MET_pheno_mean_fw$environment) %>%
  select(environment, latitude, longitude) %>%
  gather(measure, value, -environment) %>% 
  group_by(measure) %>% 
  summarize_at(vars(value), funs(min, max), na.rm = T) %>% 
  mutate(range = max - min)


### Phenotypic Stability

# Transformation

# Extract the unique stability coefficients
# Then add the breeding program information
S2_MET_pheno_fw1 <- S2_MET_pheno_mean_fw %>%
  left_join(., subset(entry_list, Class == "S2TP", c(Line, Program)), 
            by = c("line_name" = "Line")) %>%
  dplyr::rename(program = Program)

## Log transform the non-linear stability estimates
S2_MET_pheno_fw_trans <- S2_MET_pheno_fw1 %>%
  group_by(trait) %>% 
  mutate(log_delta = log(delta)) %>%
  # Tidy
  gather(term, estimate, b, delta, log_delta) %>% 
  filter(term != "delta")

# Remove potential outliers (visually)
S2_MET_pheno_fw_trans_filter <- S2_MET_pheno_fw_trans %>%
  # Remove 07MT-10 from grain protein linear stability
  filter(!(trait == "GrainProtein" & line_name == "07MT-10" & term == "b")) %>%
  # Remove 07AB-84 and 08AB-08 from heading date non-linear stability
  filter(!(trait == "HeadingDate" & term == "log_delta" & estimate > 3)) %>%
  # Remove 07MT-10 from TestWeight
  filter(!(trait == "TestWeight" & line_name == "07MT-10"))
  
S2_MET_pheno_fw_use <- S2_MET_pheno_fw_trans_filter %>%
  filter(trait %in% c("GrainYield", "HeadingDate", "PlantHeight"))

#### Summary

# Ranges and variances

# What is the range of linear stability per trait?
S2_MET_pheno_fw_use %>%
  group_by(trait) %>% 
  summarize_at(vars(h), funs(min, max)) %>% 
  mutate(range = max - min)

# Extract the distinct combinations of line, trait, and stability
S2_MET_pheno_fw_uniq <- S2_MET_pheno_fw_use %>%
  distinct(trait, line_name, g, program, term, estimate)


# What is the quartile coefficient of dispersion for each measure of stability?
qcd <- function(x) {
  quants <- quantile(x = x, probs = c(0.50, 0.75), na.rm = TRUE)
  as.numeric(diff(quants) / sum(quants))
}

pheno_mean_fw_qcd <- S2_MET_pheno_fw_uniq %>% 
  group_by(trait, term) %>%
  summarize(qcd = qcd(estimate))

# Plot
(g_qcd <- pheno_mean_fw_qcd %>% 
    ggplot(aes(x = trait, y = qcd, shape = term, color = term)) + 
    geom_point(size = 3) )


# Plot the typical FW regression (i.e. using the phenotypic means)

## Plot the stability terms individually, then combine
# First define a common list of ggplot modifiers
g_mod <- list(
  geom_density(aes(fill = "blue")),
  xlab("Estimate"),
  theme_bw() +
  theme(axis.title.y = element_blank()) )

# Just plot linear stability
g_pheno_fw_dens_b <- S2_MET_pheno_fw_uniq %>%
  filter(term == "b") %>%
  ggplot(aes(x = estimate)) + 
  labs(title = expression("Linear Stability"~(b[i]))) +
  facet_wrap( ~ trait, ncol = 1) +
  scale_fill_brewer(palette = "Set2", guide = FALSE) +
  g_mod


# Just plot non-linear stability
g_pheno_fw_dens_delta <- S2_MET_pheno_fw_uniq %>%
  filter(term == "log_delta") %>%
  ggplot(aes(x = estimate)) + 
  labs(title = expression("Non-Linear Stability"~(ln(delta[i]^2)))) +
  scale_fill_brewer(palette = "Set2", guide = FALSE) +
  facet_wrap( ~ trait, ncol = 1, scales = "free_x") +
  g_mod

# Add the plots together
g_fw_dist <- g_pheno_fw_dens_b + g_pheno_fw_dens_delta

save_file <- file.path(fig_dir, "stability_estimate_distriubtions.jpg")
ggsave(filename = save_file, plot = g_fw_dist, height = 10, width = 5, dpi = 500)




# Plot the stability estimates as lines against g_i vs h_j

# Add the program information to the results
S2_MET_pheno_mean_fw_toplot <- S2_MET_pheno_fw_use %>% 
  spread(term, estimate) %>% 
  select(trait, environment, line_name, value, g, h, b) %>% 
  left_join(., subset(entry_list, Class == "S2TP", c(Line, Program)), by = c("line_name" = "Line")) %>%
  dplyr::rename(program = Program)

# Plot the normal FW analysis plot (genotype mean against environmental effects)
g_pheno_fw_b <- S2_MET_pheno_mean_fw_toplot %>%
  ggplot(aes(x = h, y = value, col = program, group = line_name)) + 
  geom_point() + 
  geom_abline(aes(slope = b, intercept = g, col = program), alpha = 0.15) +
  facet_wrap(~ trait, ncol = 1, scales = "free") +
  scale_color_discrete(drop = FALSE, guide = guide_legend(title = "Program", ncol = 1)) +
  ylab("Genotype Mean") +
  xlab("Environmental Effect") +
  labs(title = "Linear Phenotypic Stability") +
  # scale_color_brewer(palette = "Set2") +
  theme_bw() +
  theme(legend.position = "right",
        legend.background = element_rect(fill = "grey85", color = NA))

save_file <- file.path(fig_dir, "pheno_fw_b.jpg")
ggsave(filename = save_file, plot = g_pheno_fw_b, width = 5, height = 12)


## Subset genotypes to highlight different responses
pheno_fw_example <- S2_MET_pheno_mean_fw_toplot %>%
  distinct(trait, line_name, g, b) %>%
  group_by(trait) %>% 
  filter(b == min(b) | b == max(b) | abs(1 - b) == min(abs(1 - b)))

g_pheno_fw_b_example <- S2_MET_pheno_mean_fw_toplot %>%
  select(trait, environment, line_name, value, h) %>% 
  ggplot(aes(x = h, y = value)) + 
  geom_point() + 
  geom_abline(data = pheno_fw_example, aes(intercept = g, slope = b, col = line_name), lwd = 1) + 
  scale_color_discrete(guide = guide_legend(title = "Line Name", ncol = 1)) +
  facet_wrap(~ trait, ncol = 1, scales = "free") +
  ylab("Genotype Mean") +
  xlab("Environment Effect") +
  labs(title = "Example Genotype Responses") +
  theme_bw() +
  theme(legend.position = "right",
        legend.background = element_rect(fill = "grey85", color = NA),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10))

save_file <- file.path(fig_dir, "pheno_fw_b_example.jpg")
ggsave(filename = save_file, plot = g_pheno_fw_b_example, width = 5, height = 8)


# Combine plots
g_pheno_fw_figure <- plot_grid(g_pheno_fw_b, g_pheno_fw_b_example, labels = c("A", "B"))

save_file <- file.path(fig_dir, "pheno_fw_b_combined.jpg")
ggsave(filename = save_file, plot = g_pheno_fw_figure, width = 10, height = 8)



# Count the number of cross-overs per genotype

# Group by trait
S2_MET_pheno_mean_fw_nCO <- S2_MET_pheno_mean_fw_toplot %>% 
  group_by(trait) %>% 
  do({
    df <- .
    
    # What is the range in environemntal values?
    h_range <- range(df$h)
    
    # Get the distinct values of b and g for all lines
    df1 <- distinct(df, line_name, g, b)
    
    # List of genotypes
    line_names <- unique(df1$line_name)
    
    # Create a df
    cross_over_count <- expand.grid(geno1 = line_names, geno2 = line_names) %>%
      filter(geno1 != geno2)
    
    # Iterate over the data.frame
    cross_over_count1 <- apply(X = cross_over_count, MARGIN = 1, FUN = function(i) {
        
      # Subset the information
      df1_sub <- subset(df1, line_name %in% c(i), c(b, g))
      
      # Create a square matrix
      A <- cbind(c(1, 1), -df1_sub$b)
      b <- df1_sub$g
      
      # Solve
      x <- solve(A, b)
      
      # Get the x value
      x1 <- x[2]
      
      # Is the value within the environment range?
      between(x1, h_range[1], h_range[2]) })
    
    # Combine the dfs
    cbind(cross_over_count, cross_over_count1) })


# Find the average proportion of crossovers per trait
S2_MET_pheno_mean_fw_nCO %>% 
  group_by(trait, geno1) %>% 
  summarize(nCO = sum(cross_over_count1), 
            meanCO = mean(cross_over_count1)) %>% 
  summarize(mean = mean(meanCO))

# Visualize
S2_MET_pheno_mean_fw_nCO %>% 
  group_by(trait, geno1) %>% 
  summarize(nCO = sum(cross_over_count1), 
            meanCO = mean(cross_over_count1)) %>% 
  ggplot(aes(x = trait, y = meanCO)) + 
  geom_boxplot() + 
  theme_bw()


## Plot the relationship between the genotypic effect and the sensitivity

## First fit linear models to evaluate the relationship

# Calculate the correlation between genotype mean and the stability estimate
stability_mean_corr <- S2_MET_pheno_fw_use %>%
  group_by(trait, term) %>% 
  # Bootstrap the correlation
  do(boot = {
    df <- .
    # Perform bootstrap
    boot_out <- df %>% 
      bootstrap(n = 1000) %>% 
      pull(strap) %>%
      map_dbl(~as.data.frame(.) %>% ungroup() %>% select(g, estimate) %>% cor(.) %>% .[1,2])
    
    # Summarize and export
    data.frame(cor = boot_out) %>% 
      summarize(boot_mean = mean(cor), se = sd(cor), 
                ci_lower = quantile(cor, 0.05 / 2), ci_upper = quantile(cor, 1 - (0.05 / 2))) }) %>%
  # Join with original data
  left_join(S2_MET_pheno_fw_use, .) %>%
  # Calculate the base correlation and add a character annotation
  group_by(trait, term) %>% 
  mutate(corr = cor(g, estimate, use = "complete.obs"),
         corr = str_c("r = ", round(corr, 3))) %>% 
  unnest() %>%
  ungroup()


# Fit linear models
stability_mean_lm <- stability_mean_corr %>% 
  group_by(trait, term) %>%
  do(fit = lm(g ~ estimate, data = .))

# Tidy up
stability_mean_lm_tidy <- stability_mean_lm %>%
  ungroup() %>%
  mutate(tidy_fit = map(fit, tidy)) %>%
  unnest(tidy_fit) %>% 
  filter(term1 == "estimate")


## Looks like significant regression slopes for the relationship between genotype 
## mean and $b$ for GrainYield and PlantHeight, and between genotype mean and $\delta$ 
## for HeadingDate and PlantHeight.

# The relationship between genotype mean and $\delta$ is still significant after outlier removal

# Extract the data to plot
stability_mean_corr_toplot <- stability_mean_corr %>% 
  distinct(trait, line_name, program, g, term, estimate, corr)

# Create a list of plot additions
g_add <- list(geom_point(),
              geom_smooth(method = "lm", se = FALSE, col = "black"),
              geom_text(aes(x = Inf, y = -Inf, label = corr, vjust = -1, hjust = 1.2), col = "black"),
              scale_color_discrete(name = "Program"),
              ylab("Estimate"),
              xlab("Genotype Mean"),
              theme_bw(),
              theme(legend.background = element_rect(fill = "grey85", color = NA)))

# Plot just the linear stability
g_linear_stability_and_mean <- stability_mean_corr_toplot %>%
  filter(term == "b") %>%
  ggplot(aes(x = g, y = estimate, col = program, group = FALSE)) +
  facet_wrap(~ trait, scales = "free_x", ncol = 1) +
  labs(title = "Linear Stability") +
  theme(legend.position = "right") +
  g_add
  
# Plot just the non-linear stability
g_nonlinear_stability_and_mean <- stability_mean_corr_toplot %>%
  filter(term == "log_delta") %>%
  ggplot(aes(x = g, y = estimate, col = program, group = FALSE)) +
  facet_wrap(~ trait, scales = "free", ncol = 1) +
  labs(title = "Non-Linear Stability" ) +
  theme(legend.position = "right") +
  g_add
  

# Plot grid
# Remove the legend from each plot object. Also remove the y axis title from the right-hand plot
g_plotgrid <- plot_grid(
  g_linear_stability_and_mean + theme(legend.position="none", axis.title.x = element_blank()),
  g_nonlinear_stability_and_mean + theme(legend.position="none", axis.title.x = element_blank()) + ylab(""),
  align = "h",
  labels = c("A", "B"))

# Extract the legend from one of the plots and add it back in
g_plotgrid1 <- plot_grid(
  g_plotgrid, get_legend(g_linear_stability_and_mean), rel_widths = c(2, 0.3)
)

# Add a common x-axis legend.
g_plotgrid2 <- ggdraw(add_sub(g_plotgrid1, "Genotype Mean", size = 12, hjust = 0.7))

save_file <- file.path(fig_dir, "pheno_fw_stability_mean.jpg")
ggsave(filename = save_file, plot = g_plotgrid2, width = 6, height = 6)



## Repeatability of stability using resampling

# Read in the results
load(file.path(result_dir, "S2MET_pheno_fw_resampling.RData"))

S2MET_pheno_sample_fw_tomodel <- S2MET_pheno_sample_fw %>%
  mutate(log_delta = log(delta)) %>%
  select(-delta) %>%
  gather(coef, value, b:log_delta)

# Fit a fixed-effect model
reliability_fit <- S2MET_pheno_sample_fw_tomodel %>% 
  filter(!is.infinite(value)) %>% 
  # filter(iter %in% c(1:40)) %>%
  group_by(trait, p, coef) %>% 
  do(fit = lm(value ~ line_name, .))

## Summarize
reliability_summ <- reliability_fit %>% 
  ungroup() %>% 
  mutate(anova = map(fit, ~tidy(anova(.)))) %>% 
  unnest(anova) %>% 
  group_by(trait, p, coef) %>% 
  summarize(repeatability = meansq[1] / sum(meansq))

# Plot
reliability_summ %>% 
  ggplot(aes(x = p, y = repeatability)) + 
  geom_point() + 
  facet_grid(coef~ trait) + 
  theme_bw()



## Calculate the proportion of stability variance due to genomewide markers
# First calculate relationship matrix
K <- A.mat(X = S2TP_imputed_multi_genos_mat, min.MAF = 0, max.missing = 1)

# Gather and tidy
S2_MET_pheno_fw_coef_tomodel <- S2_MET_pheno_fw_use %>%
  distinct(trait, line_name, g, term, estimate) %>%
  spread(term, estimate) %>% 
  gather(coef, estimate, g, b, log_delta)

# Now calculate heritability
SNP_herit <- S2_MET_pheno_fw_coef_tomodel %>%
  group_by(trait, coef) %>%
  do({
    df <- .
    
    df1 <- df %>% mutate(geno = line_name)

    # Fit a model in which the random genetic effect due to SNPs and the random
    # genetic effect non accounted for by SNPs are modeled
    fit2 <- relmatLmer(estimate ~ (1|line_name) + (1|geno), df1,
                       relmat = list(line_name = K))
    
    # Calculate the proportion of total phenotypic variation due to these sources
    var_prop <- VarProp(fit2)
    
    data.frame(group = c("SNP_prop", "geno_prop"), 
               prop = c(var_prop$prop[1], sum(var_prop$prop[1:2])),
               stringsAsFactors = FALSE)
    
  })

# Create a table
SNP_herit %>% 
  filter(coef != "g") %>%
  spread(group, prop)

