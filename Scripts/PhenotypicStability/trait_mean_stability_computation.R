## S2 MET Stability Computations
## 
## Use Finlay-Wilkinson regression to determine the sensitivity of lines in the
## S2MET project to the environmental mean and to environmental covariables
## 
## Use a mixed model to determine the genotype mean for each trait
## 

# Load packages and set directories
library(lme4)
library(broom)
library(modelr)
library(parallel)

# Repository directory
repo_dir <- getwd()
# Project and other directories
source(file.path(repo_dir, "source.R"))

n_cores <- detectCores()


## Calculate stability using only the TP

# For each list, what are the number of environments and number of unique lines?
env <- unique(S2_MET_BLUEs_tp$environment)
n_env <- length(env)
n_entries <- n_distinct(S2_MET_BLUEs_tp$line_name)



# Calculate for each trait
# Also pull out the random effect of each environment
S2_MET_pheno_mean <- S2_MET_BLUEs_tp %>%
  group_by(trait) %>%
  do(calc_gh(.))
    
# Ungroup
S2_MET_pheno_mean <- ungroup(S2_MET_pheno_mean)
    

  

# First create the fitted models
# Fit the model, then return df with or without outliers
S2_MET_fw_fitted <- S2_MET_pheno_mean %>%
  group_by(trait, line_name) %>%
  do(calc_stability(df = .))

# Extract the original "outliers" dataset
S2_MET_fw_fitted_outliers <- S2_MET_fw_fitted %>%
  filter(type == "outliers")


## Now we have to re-estimate the genotype mean and the environmental effect after
## outlier removal

# First, were there outliers removed?
outliers <- S2_MET_fw_fitted %>% 
  summarize(n_outlier = n[1] - n[2]) %>%
  summarize(tot_outliers = sum(n_outlier))

any_outliers <- any(outliers$tot_outliers > 0)

# Vector to record the number of outliers removed
outliers_removed <- list()

iter <- 1
outliers_removed[[iter]] <- outliers

# While loop
while(any_outliers) {
  
  # Extract the data from the outlier-removed model
  S2_MET_pheno_mean_fw <- S2_MET_fw_fitted %>% 
    filter(type == "no_outliers") %>% 
    select(trait, line_name, data) %>% 
    unnest() %>%
    select(-line_name1, -trait1)
  
  # Refit the first model
  S2_MET_pheno_mean_refit <- S2_MET_pheno_mean_fw %>%
    group_by(trait) %>%
    do({
      # Extract the data.frame
      df <- .
      
      calc_gh(df) })
  
  # Refit the second model
  S2_MET_fw_refitted <- S2_MET_pheno_mean_refit %>%
    group_by(trait, line_name) %>%
    do({
      # Grab the data
      df <- .
      
      calc_stability(df = df) })
  
  # Were there outliers removed?
  outliers <- S2_MET_fw_refitted %>% 
    summarize(n_outlier = n[1] - n[2]) %>%
    summarize(tot_outliers = sum(n_outlier))
  
  any_outliers <- any(outliers$tot_outliers > 0)
  
  # Add the outliers to the list
  iter <- iter + 1
  outliers_removed[[iter]] <- outliers
  
  S2_MET_fw_fitted <- S2_MET_fw_refitted
  
}


   
# Ungroup
fw_fitted <- ungroup(S2_MET_fw_fitted) %>%
  # Replace the "no-outlier" group with the original data
  filter(type == "no_outliers") %>% 
  bind_rows(., S2_MET_fw_fitted_outliers) %>%
  arrange(trait, line_name)

# Extract the data from the outlier-removed model
pheno_mean_fw <- S2_MET_fw_fitted %>% 
  filter(type == "no_outliers") %>% 
  select(trait, line_name, data) %>% 
  unnest() %>%
  select(-line_name1, -trait1) %>%
  ungroup()

    














## Calculate stability using the TP and the VP

# For each list, what are the number of environments and number of unique lines?
env <- unique(S2_MET_BLUEs_tpvp$environment)
# DON'T tntersect with the environments for the TP - we want different environments
n_env <- length(env)
n_entries <- n_distinct(S2_MET_BLUEs_tpvp$line_name)



# Calculate for each trait
# Also pull out the random effect of each environment
S2_MET_pheno_mean <- S2_MET_BLUEs_tpvp %>%
  group_by(trait) %>%
  do(calc_gh(.))

# Ungroup
S2_MET_pheno_mean <- ungroup(S2_MET_pheno_mean)




# First create the fitted models
# Fit the model, then return df with or without outliers
S2_MET_fw_fitted <- S2_MET_pheno_mean %>%
  group_by(trait, line_name) %>%
  do(calc_stability(df = .))

# Extract the original "outliers" dataset
S2_MET_fw_fitted_outliers <- S2_MET_fw_fitted %>%
  filter(type == "outliers")


## Now we have to re-estimate the genotype mean and the environmental effect after
## outlier removal

# First, were there outliers removed?
outliers <- S2_MET_fw_fitted %>% 
  summarize(n_outlier = n[1] - n[2]) %>%
  summarize(tot_outliers = sum(n_outlier))

any_outliers <- any(outliers$tot_outliers > 0)

# Vector to record the number of outliers removed
outliers_removed <- list()

iter <- 1
outliers_removed[[iter]] <- outliers

# While loop
while(any_outliers) {
  
  # Extract the data from the outlier-removed model
  S2_MET_pheno_mean_fw <- S2_MET_fw_fitted %>% 
    filter(type == "no_outliers") %>% 
    select(trait, line_name, data) %>% 
    unnest() %>%
    select(-line_name1, -trait1)
  
  # Refit the first model
  S2_MET_pheno_mean_refit <- S2_MET_pheno_mean_fw %>%
    group_by(trait) %>%
    do({
      # Extract the data.frame
      df <- .
      
      calc_gh(df) })
  
  # Refit the second model
  S2_MET_fw_refitted <- S2_MET_pheno_mean_refit %>%
    group_by(trait, line_name) %>%
    do({
      # Grab the data
      df <- .
      
      calc_stability(df = df) })
  
  # Were there outliers removed?
  outliers <- S2_MET_fw_refitted %>% 
    summarize(n_outlier = n[1] - n[2]) %>%
    summarize(tot_outliers = sum(n_outlier))
  
  any_outliers <- any(outliers$tot_outliers > 0)
  
  # Add the outliers to the list
  iter <- iter + 1
  outliers_removed[[iter]] <- outliers
  
  S2_MET_fw_fitted <- S2_MET_fw_refitted
  
}



# Ungroup
fw_fitted_tpvp <- ungroup(S2_MET_fw_fitted) %>%
  # Replace the "no-outlier" group with the original data
  filter(type == "no_outliers") %>% 
  bind_rows(., S2_MET_fw_fitted_outliers) %>%
  arrange(trait, line_name)

# Extract the data from the outlier-removed model
pheno_mean_fw_tpvp <- S2_MET_fw_fitted %>% 
  filter(type == "no_outliers") %>% 
  select(trait, line_name, data) %>% 
  unnest() %>%
  select(-line_name1, -trait1) %>%
  ungroup()



# Save this
save_file <- file.path(result_dir, "pheno_mean_fw_results.RData")
save("fw_fitted", "pheno_mean_fw", "fw_fitted_tpvp", "pheno_mean_fw_tpvp", file = save_file)













## What are the correlations between the environmental mean calculating using the TP versus using the TP and VP?
env_mean <- distinct(pheno_mean_fw, trait, environment, h)
env_mean_tpvp <- distinct(pheno_mean_fw_tpvp, trait, environment, h)

## Inner join
env_mean_join <- inner_join(env_mean, env_mean_tpvp, by = c("trait", "environment"))

# Correlation
(env_mean_cor <- env_mean_join %>% 
  group_by(trait) %>% 
  summarize(env_mean_cor = cor(h.x, h.y)) %>%
  mutate(annotation = str_c("r = ", formatC(x = env_mean_cor, digits = 3, format = "f"))) )

# Plot
g_env_mean_corr <- env_mean_join %>% 
  ggplot(aes(x = h.x, y = h.y)) +
  geom_point(size = 1) + 
  geom_text(data = env_mean_cor, aes(x = Inf, y = -Inf, label =  annotation), hjust = 1.5, vjust = -2, size = 2) + 
  facet_wrap(~trait, ncol = 1, scales = "free") +
  xlab("Environmental Mean (n = 183)") +
  ylab("Environmental Mean (n = 183 + 50)") +
  theme_pnas()

# Save this
ggsave(filename = "env_mean_corr.jpg", plot = g_env_mean_corr, path = fig_dir,
       height = 12, width = 8.7, units = "cm", dpi = 1000)



# trait       env_mean_cor
# 1 GrainYield         1.000
# 2 HeadingDate        0.997
# 3 PlantHeight        0.999





## How robust are our estimates of stability?
# Use a resampling approach where 20, 40, 60, or 80% of the environments are used, then estimate the stability. 
# Compare with the original estimation


# Vector of proportion of environments
p_env <- seq(0.2, 0.8, by = 0.2) %>%
  set_names(., str_c("p_env_", .))
# Number of resample iterations
n_iter <- 250



### Just the TP
# Extract data to model
pheno_tomodel <- pheno_mean_fw %>% 
  distinct(environment, line_name, trait, value, std_error, h)



# Iterate over the proportion of environments to sample
# Create bootstrapping replicates
pheno_samples <- p_env %>% 
  map(function(p) {
    # Create samples grouped by trait
    pheno_tomodel %>% 
      group_by(trait) %>% 
      do({
        trait_df <- .
        
        # Name of environments
        envs <- distinct(trait_df, environment)
        
        # Sample environments n_iter times
        sample_envs <- rerun(n_iter, sample_frac(tbl = envs, size = p))
        
        # Map over the list and subset the trait data
        sample_trait_df <- sample_envs %>% 
          map(~left_join(., trait_df, by = "environment"))
        
        # Create and return a data.frame
        data_frame(iter = seq(n_iter), data = sample_trait_df) }) %>% ungroup()
    
  }) %>% list(., names(.)) %>% pmap_df(~mutate(.x, p = .y))


# Iterate over samples and calculate the stability coefficients
pheno_samples_fit_out <- pheno_samples$data %>%
  map(~{
    df <- .
    # Group by line name and fit the models
    df %>%
      group_by(line_name) %>%
      do({
        # Fit the linear model 
        fit <- lm(value ~ h, .)
        data.frame(b = coef(fit)[2], delta = mean(resid(fit)^2), row.names = NULL, stringsAsFactors = FALSE)
      }) %>% ungroup()
  })

## Add the results to the samples df
pheno_samples_fw <- pheno_samples %>%
  mutate(p = parse_number(p),
         out = pheno_samples_fit_out) %>%
  select(-data) %>%
  unnest()
  
  
  

###  TP and VP
# Extract data to model
pheno_tomodel <- pheno_mean_fw_tpvp %>% 
  distinct(environment, line_name, trait, value, std_error, h)

# Iterate over the proportion of environments to sample
# Create bootstrapping replicates
pheno_samples <- p_env %>% 
  map(function(p) {
    # Create samples grouped by trait
    pheno_tomodel %>% 
      filter(line_name %in% tp) %>%
      group_by(trait) %>% 
      do({
        trait_df <- .
        
        # Name of environments
        envs <- distinct(trait_df, environment)
        
        # Sample environments n_iter times
        sample_envs <- rerun(n_iter, sample_frac(tbl = envs, size = p))
        
        # Map over the list and subset the trait data
        sample_trait_df <- sample_envs %>% 
          map(~left_join(., trait_df, by = "environment"))
        
        # Create and return a data.frame
        data_frame(iter = seq(n_iter), data = sample_trait_df) }) %>% ungroup()
    
  }) %>% list(., names(.)) %>% pmap_df(~mutate(.x, p = .y))


# Iterate over samples and calculate the stability coefficients
pheno_samples_fit_out <- pheno_samples$data %>%
  map(~{
    df <- .
    # Group by line name and fit the models
    df %>%
      group_by(line_name) %>%
      do({
        # Fit the linear model 
        fit <- lm(value ~ h, .)
        data.frame(b = coef(fit)[2], delta = mean(resid(fit)^2), row.names = NULL, stringsAsFactors = FALSE)
      }) %>% ungroup()
  })

## Add the results to the samples df
pheno_samples_fw_tpvp <- pheno_samples %>%
  mutate(p = parse_number(p),
         out = pheno_samples_fit_out) %>%
  select(-data) %>%
  unnest()

  

# Save the results
save_file <- file.path(result_dir, "pheno_fw_resampling.RData")
save("pheno_samples_fw","pheno_samples_fw_tpvp", file = save_file)

