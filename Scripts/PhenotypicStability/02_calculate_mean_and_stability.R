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




# Sample environments in the extremes -------------------------------------


# Only do this for traits with enough environments (GY, HD, PH)

# Run the following sampling strategies:
# 1. Upper 33%
# 2. Lower 33%
# 3. Middle 33%
# 4. Upper 33/2% and lower 33/2%
# 

# These percents are percents of the RANGE of environmental means, not the quantiles
extreme_env_sampling_probs <- list(
  upper = c(0.67, 1),
  lower = c(0, 0.33),
  middle = c(0.5 - (0.33/2), 0.5 + (0.33/2)),
  tails = c((0.33/2), 1 - (0.33/2))
)

# Number of samples for each
nRandomRep <- 10
nEnvSample <- 4


sample_extreme_envs_per_trait <- crossing(quantile = names(extreme_env_sampling_probs), rep = seq_len(nRandomRep),
                                          trait = unique(genotype_environment_mean$trait)) %>%
  filter(trait %in% c("GrainYield", "HeadingDate", "PlantHeight")) %>%
  group_by(quantile, trait, rep) %>%
  do(data = {
    row <- .
    # Get the probs to use
    probs <- extreme_env_sampling_probs[[row$quantile]]
    env_means <- subset(genotype_environment_mean, trait == row$trait)
    # Calculate the range in means to sample
    env_means_range <- diff(range(env_means$h))
    quants <- min(env_means$h) + (probs * env_means_range)
    
    if (row$quantile == "tails") {
      sample_n(tbl = filter(env_means, !between(h, quants[1], quants[2]))["environment"], size = nEnvSample)
    } else {
      sample_n(tbl = filter(env_means, between(h, quants[1], quants[2]))["environment"], size = nEnvSample)
    }
    
  }) %>% ungroup()


pheno_extreme_env_samples <- sample_extreme_envs_per_trait %>%
  mutate(data = map2(trait, data, ~subset(pheno_dat2, trait == .x & environment %in% .y$environment)))

# Calculate stability for each sample
genotype_stability_extreme_samples <- pheno_extreme_env_samples %>%
  group_by(trait, quantile, rep) %>%
  do({
    row <- .
    fwr(data = row$data[[1]], pheno = "value", gen = "line_name", env = "environment", 
        covariate = "h", weights = "smith.weights",
        model.type = "joint", genotype.random = TRUE, stability.random = TRUE)
  }) %>% ungroup()



## Analyze

# Calculate the mean absolute deviations from the stability estimates derived
# from environmental samples versus the stability estimated using all data

genotype_stability_extreme_combined <- genotype_stability_extreme_samples %>%
  mutate(stability_sample = map(regression, ~select(.x, line_name, var_d_sample = var_d, g_sample = g, b_sample = b))) %>%
  select(-regression, -ranova) %>%
  left_join(., select(genotype_stability, -ranova)) %>%
  rename(stability_allData = regression) %>%
  mutate(regression = map2(stability_sample, stability_allData, merge))

stability_extreme_sample_deviations <- genotype_stability_extreme_combined %>%
  select(trait, quantile, rep, regression) %>%
  unnest(regression) %>%
  mutate(sd_d_sample = sqrt(var_d_sample), sd_d = sqrt(var_d),
         g_sample_dev = abs(g_sample - g),
         b_sample_dev = abs(b_sample - b),
         d_sample_dev = abs(sd_d_sample - sd_d))

stability_extreme_sample_mae <- stability_extreme_sample_deviations %>%
  group_by(trait, quantile, rep) %>%
  summarize_at(vars(contains("dev")), list(mae = mean), na.rm = TRUE) %>%
  ungroup()

# Plot
stability_extreme_sample_mae %>%
  gather(parameter, sample_mae, contains("mae")) %>%
  mutate(parameter = str_remove(parameter, "_sample_dev_mae")) %>%
  filter(parameter != "g") %>%
  ggplot(aes(x = quantile, y = sample_mae)) +
  geom_boxplot() +
  facet_wrap(~ trait + parameter, ncol = 2, scales = "free_y")





# Save
save("stability_data_use", "genotype_stability_herit", "genotype_stability_cor", "genotype_environment_mean", 
     "genotype_stability", "genotype_stability_samples", "stability_sample_deviations",
     "stability_extreme_sample_deviations", "genotype_stability_extreme_samples",
     file = file.path(result_dir, "stability_estimates.RData"))

