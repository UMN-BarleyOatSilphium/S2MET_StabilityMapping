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
# Calculate for each trait
# Also pull out the random effect of each environment
S2_MET_pheno_mean <- S2_MET_BLUEs_use %>%
  group_by(trait) %>%
  do({
    # Extract the data.frame
    df <- .
    
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
    df %>% 
      left_join(., geno_blues) %>% 
      left_join(., env_blups)
    
  })
    
# Ungroup
S2_MET_pheno_mean <- ungroup(S2_MET_pheno_mean)
    


## Finlay-Wilkinson Regression
# Use the mean of each genotype in each environment and the random effect of 
# each environment to calculate the random regression coefficient of genotype
# on environment

# Group by trait and genotype
S2_MET_pheno_mean_fw <- S2_MET_pheno_mean %>%
  group_by(trait, line_name) %>%
  do({
    # Grab the data
    df <- .
    
    # Set the control
    # control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore", check.nlev.gtr.1 = "ignore")
    
    # Fit the model
    fit <- lm(value ~ h, data = df)
    # fit <- lmer(value ~ (h|line_name), data = df, control = control)
    
    # Extract the stability measures
    df %>%
      mutate(b = coef(fit)[2],
             b_std_error = subset(tidy(fit), term == "h", std.error, drop = TRUE), # Regression coefficient
             delta = mean(resid(fit)^2))
    
    # df %>% 
    #   mutate(b = ranef(fit)$line_name$h, delta = mean(resid(fit)^2))
    
  })
   
# Ungroup
S2_MET_pheno_mean_fw <- ungroup(S2_MET_pheno_mean_fw)
    
# Save this
save_file <- file.path(result_dir, "S2MET_pheno_mean_fw_results.RData")
save("S2_MET_pheno_mean_fw", file = save_file)


### Fit the environmental covariable version of the Finlay-Wilkinson Regression model

# Fit the model for the one-year and multi-year ECs
geno_stab_coef_oneyear_ec <- left_join(S2_MET_BLUEs_to_model, rename(one_year_env_df, h_ec = value), by = "environment") %>%
  group_by(trait, line_name, variable) %>%
  do({
    df <- .
    # Fit the model
    fit <- lm(value ~ h_ec, data = df)
    # Extract the coefficient table
    tidy_fit <- tidy(fit)
    data.frame(g = coef(fit)[1],
               b = coef(fit)[2], 
               b_std_error = subset(tidy_fit, term == "h_ec", std.error, drop = TRUE), 
               delta = mean(resid(fit)^2), row.names = NULL) 
    }) %>%
  gather(stability_term, estimate, -trait:-g, -b_std_error) %>%
  ungroup()

geno_stab_coef_multiyear_ec <- left_join(S2_MET_BLUEs_to_model, rename(multi_year_env_df, h_ec = value), by = "environment") %>% 
  group_by(trait, line_name, variable) %>%
  do({
    df <- .
    # Fit the model
    fit <- lm(value ~ h_ec, data = df)
    # Extract the coefficient table
    tidy_fit <- tidy(fit)
    data.frame(g = coef(fit)[1],
               b = coef(fit)[2], 
               b_std_error = subset(tidy_fit, term == "h_ec", std.error, drop = TRUE), 
               delta = mean(resid(fit)^2), row.names = NULL) 
  }) %>%
  gather(stability_term, estimate, -trait:-g, -b_std_error) %>%
  ungroup()

# Combine with the above FW results for genotype means in environments
S2_MET_ec_oneyear_fw <- S2_MET_pheno_fw %>% 
  select(., trait, environment, line_name, value, program) %>% 
  left_join(., geno_stab_coef_oneyear_ec, by = c("trait", "line_name")) %>%
  left_join(., rename(one_year_env_df, h_ec = value), by = c("environment", "variable"))

S2_MET_ec_multiyear_fw <- S2_MET_pheno_fw %>% 
  select(., trait, environment, line_name, value, program) %>% 
  left_join(., geno_stab_coef_multiyear_ec, by = c("trait", "line_name")) %>%
  left_join(., rename(multi_year_env_df, h_ec = value), by = c("environment", "variable"))

# Save  
save_file <- file.path(result_dir, "S2MET_ec_fw_regression_results.RData")
save("S2_MET_ec_oneyear_fw", "S2_MET_ec_multiyear_fw", file = save_file)






# ### Archived code for looking at both the TP and VP
# 
# ## Two populations - just the TP or the TP and the VP
# pop_all <- c(tp, vp)
# pop_tp <- tp
# 
# # Create a list of data.frames depending on the population
# S2_MET_BLUEs_to_model <- list(pop_all, pop_tp) %>%
#   map(., ~filter(S2_MET_BLUEs_use, line_name %in% .) %>%
#         droplevels() ) %>%
#   set_names(c("pop_all", "pop_tp"))
# 
# # For each list, what are the number of environments and number of unique lines?
# env <- map(S2_MET_BLUEs_to_model, ~unique(.$environment))
# n_env <- map_dbl(env, length)
# n_entries <- map_dbl(S2_MET_BLUEs_to_model, ~n_distinct(.$line_name))
# 
# 
# ## Finlay-Wilkinson Regression
# 
# # Group by trait and calculate the environmental effects
# step_one_fit <- S2_MET_BLUEs_to_model %>%
#   map(., ~group_by(., trait) %>%
#         do(step1_fit = FW(y = .$value, VAR = .$line_name, ENV = .$environment, method = "OLS")) %>%
#         ungroup() )
# 
# # Extract the environmental coefficients
# env_h_coef <- step_one_fit %>% 
#   map(., ~mutate(., h = map(step1_fit, "h") %>%
#                   map(as.data.frame) %>% 
#                   map(rownames_to_column, "environment") %>% 
#                   map(rename, h = V1)) %>% 
#         unnest(h) )
# 
# # Combine this information with the BLUEs
# # Regress each genotype and obtain the stability coefficient
# gen_stab_coef <- list(S2_MET_BLUEs_to_model, env_h_coef) %>%
#   pmap(~left_join(.x, .y, by = c("environment", "trait"))) %>%
#   map(~group_by(., trait, line_name) %>%
#         do({
#           df <- .
#           # Fit the model
#           fit <- lm(value ~ h, data = df)
#           # Extract the different stability estimates
#           tidy_fit <- tidy(fit)
#           data.frame(g = coef(fit)[1],
#                      b = coef(fit)[2], 
#                      b_std_error = subset(tidy_fit, term == "h", std.error, drop = TRUE), # Regression coefficient
#                      delta = mean(resid(fit)^2), row.names = NULL) }) %>% # delta = MSE
#         gather(stability_term, estimate, -trait:-g, -b_std_error) %>%
#         ungroup() )
# 
# # Assign entry subpopulation
# S2_MET_BLUEs_to_model_program <- S2_MET_BLUEs_to_model %>% 
#   map(~left_join(., select(entry_list, line_name = Line, program = Program), by = "line_name"))
# 
# # Levels of breeding program
# program_levels <- S2_MET_BLUEs_to_model_program %>% 
#   map(~unique(.$program))
# 
# # Add the regression coefficients to the BLUEs
# S2_MET_pheno_fw <- list(S2_MET_BLUEs_to_model_program, gen_stab_coef) %>% 
#   pmap(~left_join(.x, .y, by = c("trait", "line_name"))) %>%
#   list(., env_h_coef) %>%
#   pmap(~left_join(.x, .y, by = c("trait", "environment"))) %>% 
#   list(., program_levels) %>% 
#   pmap(~mutate(.x, program = parse_factor(program, levels = .y)))
# 
# # Save this
# save_file <- file.path(result_dir, "S2MET_pheno_fw_regression_results.RData")
# save("S2_MET_pheno_fw", file = save_file)
# 
# 
# ### Fit the environmental covariable version of the Finlay-Wilkinson Regression model
# 
# # Fit the model for the one-year and multi-year ECs
# geno_stab_coef_oneyear_ec <- S2_MET_BLUEs_to_model %>% 
#   map(., ~left_join(., rename(one_year_env_df, h_ec = value), by = "environment") %>%
#         group_by(trait, line_name, variable) %>%
#         do({
#           df <- .
#           # Fit the model
#           fit <- lm(value ~ h_ec, data = df)
#           # Extract the coefficient table
#           tidy_fit <- tidy(fit)
#           data.frame(g = coef(fit)[1],
#                      b = coef(fit)[2], 
#                      b_std_error = subset(tidy_fit, term == "h_ec", std.error, drop = TRUE), 
#                      delta = mean(resid(fit)^2), row.names = NULL) }) %>%
#         gather(stability_term, estimate, -trait:-g, -b_std_error) %>%
#         ungroup() )
#     
# geno_stab_coef_multiyear_ec <- S2_MET_BLUEs_to_model %>% 
#   map(., ~left_join(., rename(multi_year_env_df, h_ec = value), by = "environment") %>%
#         group_by(trait, line_name, variable) %>%
#         do({
#           df <- .
#           # Fit the model
#           fit <- lm(value ~ h_ec, data = df)
#           # Extract the coefficient table
#           tidy_fit <- tidy(fit)
#           data.frame(g = coef(fit)[1],
#                      b = coef(fit)[2], 
#                      b_std_error = subset(tidy_fit, term == "h_ec", std.error, drop = TRUE), 
#                      delta = mean(resid(fit)^2), row.names = NULL) }) %>%
#         gather(stability_term, estimate, -trait:-g, -b_std_error) %>%
#         ungroup() )
# 
# # Combine with the above FW results for genotype means in environments
# S2_MET_ec_oneyear_fw <- S2_MET_pheno_fw %>% 
#   map(~select(., trait, environment, line_name, value, program)) %>% 
#   list(., geno_stab_coef_oneyear_ec) %>%
#   pmap(~left_join(.x, .y, by = c("trait", "line_name"))) %>%
#   map(~left_join(., rename(one_year_env_df, h_ec = value), by = c("environment", "variable")))
# 
# S2_MET_ec_multiyear_fw <- S2_MET_pheno_fw %>% 
#   map(~select(., trait, environment, line_name, value, program)) %>% 
#   list(., geno_stab_coef_multiyear_ec) %>%
#   pmap(~left_join(.x, .y, by = c("trait", "line_name"))) %>%
#   map(~left_join(., rename(multi_year_env_df, h_ec = value), by = c("environment", "variable")))
# 
# # Save  
# save_file <- file.path(result_dir, "S2MET_ec_fw_regression_results.RData")
# save("S2_MET_ec_oneyear_fw", "S2_MET_ec_multiyear_fw", file = save_file)
# 
# 
# 
# 
# 
# 
# 
