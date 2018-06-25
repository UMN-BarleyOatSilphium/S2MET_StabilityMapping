## S2 MET Stability Computations
## 
## Use Finlay-Wilkinson regression to determine the sensitivity of lines in the
## S2MET project to the environmental mean and to environmental covariables
## 
## Use a mixed model to determine the genotype mean for each trait
## 

# Load packages and set directories
library(FW)
library(lme4)
library(broom)
library(modelr)

# Repository directory
repo_dir <- getwd()
# Project and other directories
source(file.path(repo_dir, "source.R"))



# For each list, what are the number of environments and number of unique lines?
env <- unique(S2_MET_BLUEs_use$environment)
n_env <- length(env)
n_entries <- n_distinct(S2_MET_BLUEs_use$line_name)



# Calculate for each trait
# Also pull out the random effect of each environment
S2_MET_pheno_mean <- S2_MET_BLUEs_use %>%
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

    
# Save this
save_file <- file.path(result_dir, "pheno_mean_fw_results.RData")
save("fw_fitted", "pheno_mean_fw", file = save_file)


# How robust are our estimates of stability?
#   
# Use a resampling approach where 20, 40, 60, or 80% of the environments are used, then estimate the stability. Compare with the original estimation


# Extract data to model
S2_MET_pheno_tomodel <- pheno_mean_fw %>% 
  distinct(environment, line_name, trait, value, std_error, h)

# Vector of proportion of environments
p_env <- seq(0.2, 0.8, by = 0.2)
# Number of resample iterations
n_iter <- 250

# Iterate over the proportion of environments to sample
# Create bootstrapping replicates
pheno_samples <- p_env %>% 
  set_names(., .) %>%
  map(function(p) {
    # Create samples grouped by trait
    S2_MET_pheno_tomodel %>% 
      group_by(trait) %>% 
      do({
        trait_df <- .
        
        # Name of environments
        envs <- data_frame(environment = unique(trait_df$environment))
        
        # Sample environments n_iter times
        sample_envs <- rerun(n_iter, sample_frac(tbl = envs, size = p))
        
        # Map over the list and subset the trait data
        sample_trait_df <- sample_envs %>% 
          map(~left_join(., trait_df, by = "environment"))
        
        # Create and return a data.frame
        data_frame(iter = seq(n_iter), data = sample_trait_df) }) %>% ungroup()
    
  }) %>% list(., names(.)) %>% pmap_df(~mutate(.x, p = .y))

 # Iterate over samples and calculate the stability coefficients
pheno_sample_mean_fw <- pheno_samples %>%
  unnest() %>%
  group_by(trait, p, iter, line_name) %>%
  select(-environment, -trait1, -std_error) %>%
  do({
    df <- .
    
    # Fit the linear model 
    fit <- lm(value ~ h, df)
    data.frame(b = coef(fit)[2], delta = mean(resid(fit)^2))
  }) 

# Ungroup
pheno_sample_mean_fw <- ungroup(pheno_sample_mean_fw)

# Save the results
save_file <- file.path(result_dir, "pheno_fw_resampling.RData")
save("pheno_sample_mean_fw", file = save_file)

