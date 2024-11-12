## S2MET StabilityMapping
##
## Calculate the trait mean per se, overall stability, and stability for
## specific environmental covariates
##
##

# Load packages and set directories
library(modelr)
library(lme4)

# Project and other directories
source("startup.R")


# Number of replicates of completely random environmental sampling
nRandomRep <- 25
# Range of environment proportions to sample
pEnvSample <- seq(0.1, 0.9, by = 0.1)




# Load data ---------------------------------------------------------------

load(file.path(data_dir, "project_data.RData"))

# Number of environments per trait
pheno_dat %>%
  group_by(trait) %>%
  summarize(nEnv = n_distinct(environment), .groups = "drop")


# Calculate environmental means -------------------------------------------

# Convert factors
pheno_dat1 <- pheno_dat %>%
  filter(trait %in% traits) %>%
  mutate_at(vars(environment, line_name), as.factor)


## Fit the basic model with line name as random and environment as fixed.
## Because the environments are fixed, we don't really need to worry about
## specifying the R matrix correctly
## 
## Therefore, just use LMER


# Calculate environmental index --------------------------

# Iterate over traits
genotype_environment_mean <- pheno_dat1 %>%
  group_by(trait) %>%
  do({
    df <- .
    df1 <- droplevels(df) %>%
      mutate(environment = fct_contr_sum(environment))
    
    # Fit the model
    fit <- lmer(value ~ 1 + environment + (1|line_name), data = df1)
    
    calc_gh(object = fit)
    
  }) %>% ungroup()


# Calculate overall trait stability ---------------------------------------

# Regress the phenotypic value on the environmental mean
pheno_dat2 <- pheno_dat1 %>% 
  left_join(., genotype_environment_mean) %>%
  mutate(smith.weights = 1 / (std_error^2))


# Calculate linear and non-linear stability
genotype_stability <- pheno_dat2 %>%
  group_by(trait) %>%
  do(fwr(data = ., pheno = "value", gen = "line_name", env = "environment", covariate = "h", weights = "smith.weights",
         model.type = "joint", genotype.random = TRUE, stability.random = TRUE)) %>%
  ungroup()





# Analyze stability -------------------------------------------------------



## Calculate SNP heritability of stability ##

genotype_stability1 <- genotype_stability %>%
  filter(trait %in% traits) %>%
  select(-ranova) %>%
  unnest(regression) %>%
  # Deregress g and b
  mutate(g = g / g_r2, b = b / b_r2, sd_d = sqrt(var_d_dereg)) %>%
  select(trait, line_name, g, b, sd_d)

# Iterate over traits
genotype_stability_herit <- genotype_stability1 %>%
  filter(line_name %in% tp) %>% # Only look at the TP
  gather(parameter, value, g, b, sd_d) %>%
  group_by(trait, parameter) %>%
  nest() %>%
  ungroup() %>%
  mutate(out = list(NULL))


for (i in seq_len(nrow(genotype_stability_herit))) {
  dat <- droplevels(genotype_stability_herit$data[[i]])
  
  K1 <- K[levels(dat$line_name), levels(dat$line_name)]
  
  # Create the formula 
  fixed <- value ~ 1
  random <- ~ sommer::vsr(line_name, Gu = K1)
    
  # Fit using sommer
  fit <- sommer::mmer(fixed = fixed, random = random, data = dat, verbose = FALSE)
  
  # Calculate genomic heritabilities and the standard error
  h2 <- sommer::vpredict(object = fit, transform = h2 ~ V1 / (V1 + V2))
  
  # Combine
  herit <- as_tibble(h2) %>%
    rename_all(~c("h2", "std_error"))

  # Return
  genotype_stability_herit$out[[i]] <- herit
  
} # Close loop
  
genotype_stability_herit <- genotype_stability_herit %>%
  select(-data) %>%
  unnest(out)



## Calculate correlations between mean and stability ##


# Iterate over traits
genotype_stability2 <- genotype_stability1 %>%
  filter(line_name %in% tp) %>% # Only look at the TP
  group_by(trait) %>%
  nest() %>%
  ungroup() %>%
  mutate(out = list(NULL))


for (i in seq_len(nrow(genotype_stability2))) {
  dat <- droplevels(genotype_stability2$data[[i]])
  
  K1 <- K[levels(dat$line_name), levels(dat$line_name)]
  
  # Fit bi-variate models
  bivar_combn <- combn(x = names(dat)[-1], m = 2, simplify = FALSE)
  
  # Iterate
  bivar_model_out <- list()
  for (j in seq_along(bivar_combn)) {
    
    # Create the formula 
    fixed <- reformulate(termlabels = "1", response = paste0("cbind(", paste0(bivar_combn[[j]], collapse = ","), ")"))
    random <- ~ sommer::vsr(line_name, Gu = K1)
    
    # Fit using sommer
    fit <- sommer::mmer(fixed = fixed, random = random, data = dat)
    
    # Skip if empty
    if (!is_empty(fit)) {
      
      # Calculate correlation and its standard error
      corg <- sommer::vpredict(object = fit, transform = h2 ~ V2 / sqrt((V1 * V3)))
      corr <- sommer::vpredict(object = fit, transform = h2 ~ V5 / sqrt((V4 * V6)))
      # Combine
      correl <- tibble(parameter = c("corG", "corR"),
                       estimate = c(as.numeric(corg$Estimate), as.numeric(corr$Estimate)),
                       std_error = c(as.numeric(corg$SE), as.numeric(corr$SE)))
      
      
      # Extract the genetic and residual covariance matrices
      G <- fit$sigma[[1]]
      R <- fit$sigma[[2]]
      
      
      # Output
      out <- tibble( var_pair = paste0(bivar_combn[[j]], collapse = "/"),
        P = list(cov(dat[bivar_combn[[j]]])), G = list(G), R = list(R),
        correlation = list(correl) )
      
    } else {
      # Output
      out <- tibble( var_pair = paste0(bivar_combn[[j]], collapse = "/"),
                     P = list(cov(dat[bivar_combn[[j]]])), G = list(NULL), R = list(NULL),
                    correlation = list(NULL) )
      
    }
    
    bivar_model_out[[j]] <- out
    
  } # CLose loop
  
  # Return the results
  genotype_stability2$out[[i]] <- bind_rows(bivar_model_out)
  
}


# Unnest
genotype_stability_cor <- genotype_stability2 %>%
  select(-data) %>%
  unnest(out)








# Resample environments to calculate stability ----------------------------

## Random resampling

# Take a list of distinct environments per trait and randomly sample them
trait_envs <- distinct(pheno_dat2, trait, environment)

sample_envs_per_trait <- crossing(pEnv = pEnvSample, rep = seq_len(nRandomRep)) %>%
  mutate(samples = map(pEnv, ~ungroup(sample_frac(tbl = group_by(trait_envs, trait), size = .x)))) %>% 
  unnest(samples) %>%
  group_by(trait, pEnv, rep) %>% 
  nest() %>% 
  ungroup() %>%
  mutate(nEnv = map_dbl(data, nrow)) %>%
  # Remove cases where the number of envs is < 3
  filter(nEnv >= 3)

pheno_env_samples <- sample_envs_per_trait %>%
  mutate(data = map2(trait, data, ~subset(pheno_dat2, trait == .x & environment %in% .y$environment)))

# Calculate stability for each sample
genotype_stability_samples <- pheno_env_samples %>%
  group_by(trait, pEnv, nEnv, rep) %>%
  do({
    row <- .
    fwr(data = row$data[[1]], pheno = "value", gen = "line_name", env = "environment", 
        covariate = "h", weights = "smith.weights",
        model.type = "joint", genotype.random = TRUE, stability.random = TRUE)
  }) %>% ungroup()



# Rename pheno_dat2
stability_data_use <- pheno_dat2




# Analyze stability samples -----------------------------------------------


# Calculate the mean absolute deviations from the stability estimates derived
# from environmental samples versus the stability estimated using all data

genotype_stability_combined <- genotype_stability_samples %>%
  mutate(stability_sample = map(regression, ~select(.x, line_name, var_d_sample = var_d, g_sample = g, b_sample = b))) %>%
  select(-regression, -ranova) %>%
  left_join(., select(genotype_stability, -ranova)) %>%
  rename(stability_allData = regression) %>%
  mutate(regression = map2(stability_sample, stability_allData, merge))

stability_sample_deviations <- genotype_stability_combined %>%
  select(trait, pEnv, nEnv, rep, regression) %>%
  unnest(regression) %>%
  mutate(sd_d_sample = sqrt(var_d_sample), sd_d = sqrt(var_d),
         g_sample_dev = abs(g_sample - g),
         b_sample_dev = abs(b_sample - b),
         d_sample_dev = abs(sd_d_sample - sd_d))
         
stability_sample_mae <- stability_sample_deviations %>%
  group_by(trait, pEnv, nEnv, rep) %>%
  summarize_at(vars(contains("dev")), list(mae = mean), na.rm = TRUE) %>%
  ungroup()

# Plot
stability_sample_mae %>%
  gather(parameter, sample_mae, contains("mae")) %>%
  mutate(parameter = str_remove(parameter, "_sample_dev_mae")) %>%
  filter(parameter != "g") %>%
  ggplot(aes(x = pEnv, y = sample_mae)) +
  geom_boxplot(aes(group = pEnv)) +
  scale_x_continuous(breaks = pretty) +
  facet_wrap(~ trait + parameter, ncol = 4, scales = "free_y")




# Targeted sampling of environments -------------------------------------------


# Only do this for traits with enough environments (GY, HD, PH)

# Number of TP environments per trait
tp_envs <- stability_data_use %>%
  group_by(trait, environment) %>%
  summarize(nInd = n_distinct(line_name), .groups = "drop") %>%
  subset(nInd >= 150) 

tp_traits <- tp_envs %>%
  group_by(trait) %>%
  summarize(nEnv = n_distinct(environment), .groups = "drop")


# Run the following sampling strategies:
# 1. Upper 50%
# 2. Lower 50%
# 3. Middle 50%
# 4. Upper 25% and lower 25%
# 

# These percents are percents of the RANGE of environmental means, not the quantiles
extreme_env_sampling_probs <- list(
  upper = c(0.50, 1),
  lower = c(0, 0.50),
  middle = c(0.5 - (0.5/2), 0.5 + (0.5/2)),
  tails = c((0.5/2), 1 - (0.5/2))
)

# New proportions since we are only sampling 50% of the environments
pEnvSample1 <- pEnvSample * 2
pEnvSample1 <- subset(pEnvSample1, pEnvSample1 <= 1.0)
names(pEnvSample1) <- paste0("prop", pEnvSample1/2)

# Number of samples for each

sample_extreme_envs_per_trait <- crossing(quantile = names(extreme_env_sampling_probs), 
                                          prop = pEnvSample1,
                                          rep = seq_len(nRandomRep),
                                          trait = subset(tp_traits, nEnv >= 30, trait, drop = TRUE)) %>%
  group_by(quantile, prop, trait, rep) %>%
  do(data = {
    row <- .
    # Get the probs to use
    probs <- extreme_env_sampling_probs[[row$quantile]]
    env_means <- subset(genotype_environment_mean, trait == row$trait) %>%
      inner_join(., tp_envs, by = c("trait", "environment"))
    
    quants <- quantile(env_means$h, probs = probs)
    
    # Subset environments within the quantile
    env_means_quant <- env_means[between(env_means$h, quants[1], quants[2]),]
    # Randomly sample these environments
    env_means_quant[sort(sample(nrow(env_means_quant), ceiling(row$prop * nrow(env_means_quant)))),]

  }) %>% ungroup()


# No need to re-calculate environmental means because they will be the same
pheno_extreme_env_samples <- sample_extreme_envs_per_trait %>%
  mutate(data = map2(trait, data, ~subset(pheno_dat2, trait == .x & environment %in% .y$environment)))

genotype_stability_extreme_samples <- list()
for (i in seq_len(nrow(pheno_extreme_env_samples))) {
  row <- pheno_extreme_env_samples[i,]
  fwr_out <- fwr(data = row$data[[1]], pheno = "value", gen = "line_name", env = "environment", 
                 covariate = "h", weights = "smith.weights", 
                 model.type = "joint", genotype.random = TRUE, stability.random = TRUE)
  genotype_stability_extreme_samples[[i]] <- fwr_out
}

genotype_stability_extreme_samples_df <- pheno_extreme_env_samples %>%
  mutate(results = genotype_stability_extreme_samples) %>%
  select(-data) %>%
  unnest(results)


## Analyze

# Calculate the mean absolute deviations from the stability estimates derived
# from environmental samples versus the stability estimated using all data

genotype_stability_extreme_combined <- genotype_stability_extreme_samples_df %>%
  mutate(stability_sample = map(regression, ~select(.x, line_name, var_d_sample = var_d, g_sample = g, b_sample = b))) %>%
  dplyr::select(-regression, -ranova) %>%
  left_join(., dplyr::select(genotype_stability, -ranova, -vars)) %>%
  dplyr::rename(stability_allData = regression) %>%
  mutate(regression = map2(stability_sample, stability_allData, merge))

stability_extreme_sample_deviations <- genotype_stability_extreme_combined %>%
  select(trait, quantile, prop, rep, regression) %>%
  unnest(regression) %>%
  mutate(sd_d_sample = sqrt(var_d_sample), sd_d = sqrt(var_d),
         g_sample_dev = abs(g_sample - g),
         b_sample_dev = abs(b_sample - b),
         d_sample_dev = abs(sd_d_sample - sd_d))

stability_extreme_sample_mae <- stability_extreme_sample_deviations %>%
  group_by(trait, quantile, prop, rep) %>%
  summarize_at(vars(contains("dev")), list(mae = mean), na.rm = TRUE) %>%
  ungroup()

# Plot
stability_extreme_sample_mae1 <- stability_extreme_sample_mae %>%
  gather(parameter, sample_mae, contains("mae")) %>%
  mutate(parameter = str_remove(parameter, "_sample_dev_mae"))

# summarize
stability_extreme_sample_mae_summ <- stability_extreme_sample_mae1 %>%
  group_by(trait, quantile, prop, parameter) %>%
  summarize_at(vars(sample_mae), list(mean = mean, q90 = ~ quantile(., 0.9), q10 = ~ quantile(., 0.10))) %>%
  ungroup()

g_stab_extreme <- stability_extreme_sample_mae_summ %>%
  subset(parameter != "d") %>%
  ggplot(aes(x = prop, y = mean, linetype = quantile)) +
  geom_ribbon(aes(ymin = q10, ymax = q90, linetype = quantile), alpha = 0.75, fill = "grey85") +
  geom_line() +
  facet_wrap(~ trait + parameter, scales = "free_y", ncol = 2, labeller = labeller(trait = str_add_space, parameter = coef_replace1, .multi_line = FALSE)) +
  scale_x_continuous(name = "Proportion of environments sampled from the quantile", breaks = pretty) +
  scale_y_continuous(name = "Mean absolute error", breaks = pretty, limits = c(0, NA)) +
  scale_linetype_discrete(name = "Quantile of\nenvironment\nmean", labels = c("lower" = "lower 50%", "upper" = "upper 50%", "middle" = "middle 50%", "tails" = "upper/lower 25%")) +
  theme_presentation(12) +
  theme(panel.grid.major.y = element_line(), panel.spacing.y = unit(1, "line"), strip.placement = "outside")

# SAve as figure
ggsave(filename = "figureS3_extreme_environment_parameter_mae.jpg", plot = g_stab_extreme, path = fig_dir, width = 6, height = 6, dpi = 1000)
  

# Save
save("stability_data_use", "genotype_stability_herit", "genotype_stability_cor", "genotype_environment_mean", 
     "genotype_stability", "genotype_stability_samples", "stability_sample_deviations",
     "stability_extreme_sample_deviations", "genotype_stability_extreme_samples_df",
     file = file.path(result_dir, "stability_estimates.RData"))

