## S2 MET Stability Computations
## 
## Use Finlay-Wilkinson regression to determine the sensitivity of lines in the
## S2MET project to the environmental mean and to environmental covariables
## 

# Load packages and set directories

library(tidyverse)
library(readxl)
library(lme4)
library(broom)
library(stringr)
library(FW)
library(ggridges)
library(pbr)

# Project and other directories
source("C:/Users/Jeff/Google Drive/Barley Lab/Projects/S2MET_Mapping/source.R")


# Remove the environments in which the vp was only observed
S2_MET_BLUEs_use <- S2_MET_BLUEs %>%
  group_by(trait, environment) %>%
  filter(sum(line_name %in% tp) > 1) %>%
  ungroup()

# Number of environments and entries
env <- pull(distinct(S2_MET_BLUEs_use, environment))
n_env <- n_distinct(env)
n_entries <- n_distinct(S2_MET_BLUEs_use$line_name)



## Regression

### Fit the normal Finlay-Wilkinson Regression model


# Group by trait and fit the model for step 1
step_one_fit <- S2_MET_BLUEs_use %>% 
  group_by(trait) %>%
  do(step1_fit = FW(y = .$value, VAR = .$line_name, ENV = .$environment, method = "OLS")) %>%
  ungroup()

# Extract the environmental coefficients
env_h_coef <- step_one_fit %>% 
  mutate(h = map(step1_fit, pluck, "h") %>% 
           map(as.data.frame) %>% 
           map(rownames_to_column, "environment") %>% 
           map(rename, h = V1)) %>% 
  unnest(h)

# Combine this information with the BLUEs
# Regress each genotype and obtain the stability coefficient
gen_coef <- S2_MET_BLUEs_use %>% 
  left_join(., env_h_coef) %>%
  group_by(trait, line_name) %>%
  do({
    df <- .
    # Fit the model
    fit <- lm(value ~ h, data = df)
    # Extract the coefficient table
    tidy_fit <- tidy(fit)
    data.frame(g = coef(fit)[1], b = coef(fit)[2], 
               b_std_error = subset(tidy_fit, term == "h", std.error, drop = TRUE), 
               delta = sigma(fit)^2, row.names = NULL) }) %>%
  gather(stability_term, estimate, -trait:-g, -b_std_error) %>%
  ungroup()

# Color by entry subpopulation (i.e. program)
S2_MET_BLUEs_use_entries <- left_join(S2_MET_BLUEs_use, select(entry_list, Line, Program), 
                                      by = c("line_name" = "Line"))

# Levels of breeding program
program_levels <- sort(unique(S2_MET_BLUEs_use_entries$Program))

# Add the regression coefficients
S2_MET_BLUEs_fw <- S2_MET_BLUEs_use_entries %>% 
  left_join(., gen_coef, by = c("trait", "line_name")) %>% 
  left_join(., env_h_coef, by = c("trait", "environment")) %>%
  ungroup() %>%
  mutate(Program = parse_factor(Program, levels = program_levels))

# Save this
save_file <- file.path(result_dir, "S2MET_fw_regression_results.RData")
save("S2_MET_BLUEs_fw", file = save_file)


### Fit the environmental covariable version of the Finlay-Wilkinson Regression model


# Fit the model for the one-year and multi-year ECs
gen_coef_one_year <- S2_MET_BLUEs_use %>% 
  left_join(., rename(one_year_env_df, h_ec = value)) %>%
  group_by(trait, line_name, variable) %>%
  do({
    df <- .
    # Fit the model
    fit <- lm(value ~ h_ec, data = df)
    # Extract the coefficient table
    tidy_fit <- tidy(fit)
    data.frame(intercept = coef(fit)[1],
               b = coef(fit)[2], 
               b_std_error = subset(tidy_fit, term == "h_ec", std.error, drop = TRUE), 
               delta = sigma(fit)^2, row.names = NULL) }) %>%
  gather(stability_term, estimate, -trait:-intercept, -b_std_error) %>%
  ungroup()
    
gen_coef_multi_year <- S2_MET_BLUEs_use %>% 
  left_join(., rename(multi_year_env_df, h_ec = value)) %>%
  group_by(trait, line_name, variable) %>%
  do({
    df <- .
    # Fit the model
    fit <- lm(value ~ h_ec, data = df)
    # Extract the coefficient table
    tidy_fit <- tidy(fit)
    data.frame(intercept = coef(fit)[1],
               b = coef(fit)[2], 
               b_std_error = subset(tidy_fit, term == "h_ec", std.error, drop = TRUE), 
               delta = sigma(fit)^2, row.names = NULL) }) %>%
  gather(stability_term, estimate, -trait:-intercept, -b_std_error) %>%
  ungroup()

# Combine with the above FW results
S2_MET_BLUEs_one_year_fw <- S2_MET_BLUEs_fw %>% 
  select(trait, environment, line_name, value, Program, g) %>% 
  left_join(., gen_coef_one_year) %>%
  left_join(., rename(one_year_env_df, h_ec = value))

S2_MET_BLUEs_multi_year_fw <- S2_MET_BLUEs_fw %>% 
  select(trait, environment, line_name, value, Program, g) %>% 
  left_join(., gen_coef_multi_year) %>%
  left_join(., rename(multi_year_env_df, h_ec = value))

# Save  
save_file <- file.path(result_dir, "S2MET_ec_fw_regression_results.RData")
save("S2_MET_BLUEs_one_year_fw", "S2_MET_BLUEs_multi_year_fw", file = save_file)







