## S2 MET Stability Computations
## 
## Use Finlay-Wilkinson regression to determine the sensitivity of lines in the
## S2MET project to the environmental mean and to environmental covariables
## 
## Use a mixed model to determine the genotype mean for each trait
## 

# Load packages and set directories

library(tidyverse)
library(readxl)
library(lme4)
library(modelr)
library(broom)
library(stringr)
library(FW)

# Project and other directories
source("C:/Users/Jeff/Google Drive/Barley Lab/Projects/S2MET_Mapping/source.R")

tp_geno <- tp_geno_multi

# Remove the environments in which only the vp was observed
S2_MET_BLUEs_use <- S2_MET_BLUEs %>%
  filter(line_name %in% tp_geno,
         trait != "TestWeight")

# How many environments?
n_distinct(S2_MET_BLUEs_use$environment)
# Per trait?
S2_MET_BLUEs_use %>% 
  group_by(trait) %>% 
  summarize(n_env = n_distinct(environment))



# For each list, what are the number of environments and number of unique lines?
env <- unique(S2_MET_BLUEs_use$environment)
n_env <- length(env)
n_entries <- n_distinct(S2_MET_BLUEs_use$line_name)

## Genotype mean
# Function to calculate the genotype mean and the environmental effect
calc_gh <- function(df) {
  
  # Set the control
  control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
  # Extract the weights - square of the standard error of the mean
  wts <- df$std_error^2
  
  # First calculate the least-square means of each genotype
  # Fit environment and GxE as random effects
  fit1 <- lmer(value ~ -1 + line_name + (1|environment),
               data = df, control = control, weights = wts)
  
  # Extract the BLUEs
  geno_blues <- tidy(fit1) %>% 
    filter(group == "fixed") %>% 
    mutate(line_name = str_replace(term, "line_name", "")) %>% 
    select(line_name, g = estimate)
  
  # Extract the BLUPs
  env_blups <- ranef(fit1) %>% 
    as.data.frame() %>% 
    filter(grpvar == "environment") %>% 
    select(environment = grp, h = condval)
  
  # Return
  df1 <- df %>% 
    left_join(., geno_blues) %>% 
    left_join(., env_blups)
  
  return(df1)
  
}


# Calculate for each trait
# Also pull out the random effect of each environment
S2_MET_pheno_mean <- S2_MET_BLUEs_use %>%
  group_by(trait) %>%
  do({
    # Extract the data.frame
    df <- .
    
    calc_gh(df) })
    
# Ungroup
S2_MET_pheno_mean <- ungroup(S2_MET_pheno_mean)
    


## Finlay-Wilkinson Regression
# Use the mean of each genotype in each environment and the random effect of 
# each environment to calculate the random regression coefficient of genotype
# on environment


# Create a function to iteratively remove outliers based on studentized residuals
# Returns a df to be used to fit the final model
remove_outliers <- function(df, fit, cutoff = 3) {
  
  # Add residuals to the df
  df_resid <- df %>% 
    add_residuals(fit) %>%
    mutate(stand_resid = as.numeric(scale(resid, center = FALSE)))
  
  # If there are std residuals that are greater than the cutoff, remove them
  # and refit
  while (any(abs(df_resid$stand_resid) > cutoff)) {
    
    # Filter
    df_filter <- df_resid %>%
      filter(abs(stand_resid) <= cutoff)
    
    # Refit
    re_fit <- lm(value ~ h, data = df_filter)
    
    # Add the residuals and repeat
    df_resid <- df_filter %>% 
      add_residuals(re_fit) %>%
      mutate(stand_resid = as.numeric(scale(resid, center = FALSE)))
  }
  
  # Return the data.frame
  return(df_resid)
}

# Function to estimate the slope and MSE of a line
# Incorportate code to remove outliers
calc_stability <- function(df) {
  
  # # Set the control
  # control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore", check.nlev.gtr.1 = "ignore")
  # # Get the weights
  # wts <- df$std_error^2
  # fit <- lmer(value ~ (h|line_name), data = df, control = control, weights = wts)
  
  # Fit the model and return
  fit <- lm(value ~ h, data = df) 
  
  # Add residuals to the df
  # Filter for outliers
  df_filter <- remove_outliers(df = df, fit = fit, cutoff = 3)
  
  # Refit
  re_fit <- lm(value ~ h, data = df_filter) 
  
  # Return a data.frame with each data.frame, slope and MSE estimates,
  # and number of observations
  df_fit <- df %>%
    mutate(b = coef(fit)[2],
           b_std_error = subset(tidy(fit), term == "h", std.error, drop = TRUE), # Regression coefficient
           delta = mean(resid(fit)^2))
  
  df_re_fit <- df_filter %>%
    select(-contains("resid")) %>%
    mutate(b = coef(re_fit)[2],
           b_std_error = subset(tidy(re_fit), term == "h", std.error, drop = TRUE), # Regression coefficient
           delta = mean(resid(re_fit)^2))
  
  # List of data.frame
  df_list <- list(df_fit, df_re_fit)
  
  # Return data
  df1 <- data_frame(type = c("outliers", "no_outliers"), 
             data = df_list, 
             model = list(fit, re_fit),
             n = map_dbl(df_list, nrow), 
             b = map_dbl(df_list, ~unique(.$b)), 
             delta = map_dbl(df_list, ~unique(.$delta)))
  
  return(df1)
  
}

  

# First create the fitted models
# Fit the model, then return df with or without outliers
S2_MET_fw_fitted <- S2_MET_pheno_mean %>%
  group_by(trait, line_name) %>%
  do({
    # Grab the data
    df <- .
    
    calc_stability(df = df) })

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
    select(trait:z_score) %>%
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
S2_MET_fw_fitted <- ungroup(S2_MET_fw_fitted) %>%
  # Replace the "no-outlier" group with the original data
  filter(type == "no_outliers") %>% 
  bind_rows(., S2_MET_fw_fitted_outliers) %>%
  arrange(trait, line_name)

# Extract the data from the outlier-removed model
S2_MET_pheno_mean_fw <- S2_MET_fw_fitted %>% 
  filter(type == "no_outliers") %>% 
  select(trait, line_name, data) %>% 
  unnest() %>%
  select(-line_name1, -trait1)

    
# Save this
save_file <- file.path(result_dir, "S2MET_pheno_mean_fw_results.RData")
save("S2_MET_fw_fitted", "S2_MET_pheno_mean_fw", file = save_file)


# How robust are our estimates of stability?
#   
# Use a resampling approach where 20, 40, 60, or 80% of the environments are used, then estimate the stability. Compare with the original estimation


# Extract data to model
S2_MET_pheno_tomodel <- S2_MET_pheno_mean_fw %>% 
  distinct(environment, line_name, trait, value, std_error, h)

# Vector of proportion of environments
p_env <- seq(0.2, 0.8, by = 0.2)
# Number of bootstrap iterations
n_iter <- 250

# Iterate over the proportion of environments to sample
# Create bootstrapping replicates
S2_MET_pheno_samples <- p_env %>% 
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
S2MET_pheno_sample_fw <- S2_MET_pheno_samples %>%
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
S2MET_pheno_sample_fw <- ungroup(S2MET_pheno_sample_fw)

# Save the results
save_file <- file.path(result_dir, "S2MET_pheno_fw_resampling.RData")
save("S2MET_pheno_sample_fw", file = save_file)





# ### Fit the environmental covariable version of the Finlay-Wilkinson Regression model
# 
# # Fit the model for the one-year and multi-year ECs
# geno_stab_coef_oneyear_ec <- left_join(S2_MET_BLUEs_to_model, rename(one_year_env_df, h_ec = value), by = "environment") %>%
#   group_by(trait, line_name, variable) %>%
#   do({
#     df <- .
#     # Fit the model
#     fit <- lm(value ~ h_ec, data = df)
#     # Extract the coefficient table
#     tidy_fit <- tidy(fit)
#     data.frame(g = coef(fit)[1],
#                b = coef(fit)[2], 
#                b_std_error = subset(tidy_fit, term == "h_ec", std.error, drop = TRUE), 
#                delta = mean(resid(fit)^2), row.names = NULL) 
#     }) %>%
#   gather(stability_term, estimate, -trait:-g, -b_std_error) %>%
#   ungroup()
# 
# geno_stab_coef_multiyear_ec <- left_join(S2_MET_BLUEs_to_model, rename(multi_year_env_df, h_ec = value), by = "environment") %>% 
#   group_by(trait, line_name, variable) %>%
#   do({
#     df <- .
#     # Fit the model
#     fit <- lm(value ~ h_ec, data = df)
#     # Extract the coefficient table
#     tidy_fit <- tidy(fit)
#     data.frame(g = coef(fit)[1],
#                b = coef(fit)[2], 
#                b_std_error = subset(tidy_fit, term == "h_ec", std.error, drop = TRUE), 
#                delta = mean(resid(fit)^2), row.names = NULL) 
#   }) %>%
#   gather(stability_term, estimate, -trait:-g, -b_std_error) %>%
#   ungroup()
# 
# # Combine with the above FW results for genotype means in environments
# S2_MET_ec_oneyear_fw <- S2_MET_pheno_fw %>% 
#   select(., trait, environment, line_name, value, program) %>% 
#   left_join(., geno_stab_coef_oneyear_ec, by = c("trait", "line_name")) %>%
#   left_join(., rename(one_year_env_df, h_ec = value), by = c("environment", "variable"))
# 
# S2_MET_ec_multiyear_fw <- S2_MET_pheno_fw %>% 
#   select(., trait, environment, line_name, value, program) %>% 
#   left_join(., geno_stab_coef_multiyear_ec, by = c("trait", "line_name")) %>%
#   left_join(., rename(multi_year_env_df, h_ec = value), by = c("environment", "variable"))
# 
# # Save  
# save_file <- file.path(result_dir, "S2MET_ec_fw_regression_results.RData")
# save("S2_MET_ec_oneyear_fw", "S2_MET_ec_multiyear_fw", file = save_file)


