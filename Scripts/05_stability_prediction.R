## S2MET GWAS of Stability Coefficients
##
## Run genomewide predictions of the genotype mean, linear stability, and non-linear stability
## Use both within population CV and cross-population validation
## Compare multi-variate vs univariate models
## Also look at the impact of environmental sampling on prediction accuracy
## 
##
##

# Load packages and run the startup script
library(sommer)

source(here::here("startup.R"))

library(modelr) # For generating cross-validation sets

# Number of replicates for random CV
cv_reps <- 100
# Number of folds for k-fold CV
cv_folds <- 5

# Alpha level
a <- 0.10



# Load data ---------------------------------------------------------------


# Load project data
load(file.path(data_dir, "project_data.RData"))

## Load stability data
load(file.path(result_dir, "stability_estimates.RData"))

# Data.frame of stability samples
genotype_stability_df <- genotype_stability %>%
  filter(trait %in% traits) %>%
  select(-ranova) %>%
  unnest(regression) %>%
  mutate(sd_d = sqrt(var_d)) %>%
  select(trait, line_name, g, b, sd_d) %>%
  gather(parameter, value, g, b, sd_d)




# Run k-fold cross-validation ---------------------------------------------

# Create replicated training/test sets
genotype_stability_kfold_train_test <- genotype_stability_df %>%
  filter(line_name %in% tp) %>%
  group_by(trait, parameter) %>%
  do(imap_dfr(replicate(n = cv_reps, crossv_kfold(data = ., k = cv_folds), simplify = FALSE), ~mutate(.x, .rep = .y))) %>%
  ungroup()

# Run cross-val
genotype_stability_kfold_results <- genotype_stability_kfold_train_test %>%
  mutate(results = map2(train, test, ~gblup_crossv(fixed = value ~ 1, random = ~ line_name, K = K, .train = .x, .test = .y))) %>%
  select(trait, parameter, .id, .rep, results)

# Calculate accuracy
genotype_stability_kfold_acc <- genotype_stability_kfold_results %>%
  mutate(predictions = map(results, "predictions")) %>%
  select(.id, .rep, predictions) %>%
  unnest(predictions) %>%
  group_by(trait, parameter, .rep) %>%
  summarize(accuracy = cor(value, predicted.value), rmse = sqrt(mean((value - predicted.value)^2)), .groups = "drop")
  

# Plot
genotype_stability_kfold_acc %>%
  ggplot(aes(x = trait, y = accuracy, color = parameter)) +
  geom_boxplot()




# Run between-generation prediction ---------------------------------------


# Create replicated training/test sets
genotype_stability_pov_train_test <- genotype_stability_df %>%
  mutate(group = ifelse(line_name %in% tp, "train", "test")) %>%
  group_by(trait, parameter, group) %>%
  nest() %>%
  ungroup() %>%
  spread(group, data)

# Run pov
genotype_stability_pov_results <- genotype_stability_pov_train_test %>%
  mutate(results = map2(train, test, ~gblup_crossv(fixed = value ~ 1, random = ~ line_name, K = K, .train = .x, .test = .y))) %>%
  select(trait, parameter, results)

# Calculate accuracy
genotype_stability_pov_acc <- genotype_stability_pov_results %>%
  mutate(predictions = map(results, "predictions")) %>%
  unnest(predictions) %>%
  group_by(trait, parameter) %>%
  summarize(accuracy = cor(value, predicted.value), rmse = sqrt(mean((value - predicted.value)^2)), .groups = "drop")


# Plot
genotype_stability_pov_acc %>%
  ggplot(aes(x = trait, y = accuracy, color = parameter)) +
  geom_point()




# Run between-generation prediction using samples -------------------------


# The training sets will be the samples for the FP;
# the testing set will be the stability estimates of the OP using the full dataset

# Data.frame of stability samples
genotype_stability_samples_df <- genotype_stability_samples %>%
  filter(trait %in% traits) %>%
  select(-ranova, -data) %>%
  unnest(regression) %>%
  mutate(sd_d = sqrt(var_d)) %>%
  select(trait, line_name, pEnv, rep, g, b, sd_d) %>%
  gather(parameter, value, g, b, sd_d)


# Create replicated training/test sets
genotype_stability_samples_pov_train_test <- genotype_stability_samples_df %>%
  mutate(group = ifelse(line_name %in% tp, "train", "validate")) %>%
  group_by(trait, pEnv, rep, parameter, group) %>%
  nest() %>%
  ungroup() %>%
  spread(group, data) %>%
  left_join(., select(genotype_stability_pov_train_test, -train))

# Run pov
genotype_stability_samples_pov_results <- genotype_stability_samples_pov_train_test %>%
  filter(!map_lgl(genotype_stability_samples_pov_train_test$train, ~all(is.na(.x$value)))) %>%
  rowid_to_column() %>%
  group_by(trait, pEnv, rep, parameter) %>%
  do({
    row <- .
    # print(row$rowid)
    val <- gblup_crossv(fixed = value ~ 1, random = ~ line_name, K = K, .train = row$train[[1]], .test = row$validate[[1]])
    tst <- gblup_crossv(fixed = value ~ 1, random = ~ line_name, K = K, .train = row$train[[1]], .test = row$test[[1]])
    tibble(results_validation = list(val), results_test = list(tst))
  }) %>% ungroup()
  
# Calculate accuracy
genotype_stability_samples_pov_acc <- genotype_stability_samples_pov_results %>%
  mutate(predictions_validation = map(results_validation, "predictions"),
         predictions_test = map(results_test, "predictions")) %>%
  select(-contains("results")) %>%
  gather(target, predictions, contains("predictions")) %>%
  unnest(predictions) %>%
  group_by(trait, parameter, pEnv, rep, target) %>%
  summarize(accuracy = cor(value, predicted.value), rmse = sqrt(mean((value - predicted.value)^2)), .groups = "drop")



genotype_stability_samples_pov_acc_summ <- genotype_stability_samples_pov_acc %>%
  group_by(trait, parameter, pEnv, target) %>%
  summarize_at(vars(accuracy), list(mean = mean, upper = ~quantile(., 1 - a/2, na.rm = T), lower = ~quantile(., a/2, na.rm = T)), 
               na.rm = TRUE) %>%
  ungroup() %>%
  mutate(group = paste0(pEnv, target, sep  = "."))

# Plot


genotype_stability_samples_pov_acc_summ %>%
  ggplot(aes(x = pEnv, y = mean, color = target)) +
  # geom_ribbon(aes(ymin = lower, ymax = upper, fill = target), alpha = 0.2) +
  geom_line() +
  facet_grid(trait ~ parameter)




# Run between-generation prediction using extreme samples -------------------------


# The training sets will be the samples for the FP;
# the testing set will be the stability estimates of the OP using the full dataset

# Data.frame of stability samples
genotype_stability_extreme_samples_df1 <- genotype_stability_extreme_samples_df %>%
  filter(trait %in% traits) %>%
  select(-ranova) %>%
  unnest(regression) %>%
  mutate(sd_d = sqrt(var_d)) %>%
  select(trait, line_name, quantile, pEnv = prop, rep, g, b, sd_d) %>%
  gather(parameter, value, g, b, sd_d)


# Create replicated training/test sets
genotype_stability_extreme_samples_pov_train_test <- genotype_stability_extreme_samples_df1 %>%
  mutate(group = ifelse(line_name %in% tp, "train", "validate")) %>%
  group_by(trait, quantile, pEnv, rep, parameter, group) %>%
  nest() %>%
  ungroup() %>%
  spread(group, data) %>%
  left_join(., select(genotype_stability_pov_train_test, -train))

# Run pov
genotype_stability_extreme_samples_pov_results <- genotype_stability_extreme_samples_pov_train_test %>%
  filter(!map_lgl(genotype_stability_extreme_samples_pov_train_test$train, ~all(is.na(.x$value)))) %>%
  # subset(parameter == "b") %>%
  rowid_to_column() %>%
  group_by(trait, quantile, pEnv, rep, parameter) %>%
  do({
    row <- .
    # print(row$rowid)
    val <- gblup_crossv(fixed = value ~ 1, random = ~ line_name, K = K, .train = row$train[[1]], .test = row$validate[[1]])
    tst <- gblup_crossv(fixed = value ~ 1, random = ~ line_name, K = K, .train = row$train[[1]], .test = row$test[[1]])
    tibble(results_validation = list(val), results_test = list(tst))
  }) %>% ungroup()

# Calculate accuracy
genotype_stability_extreme_samples_pov_acc <- genotype_stability_extreme_samples_pov_results %>%
  mutate(predictions_validation = map(results_validation, "predictions"),
         predictions_test = map(results_test, "predictions")) %>%
  select(-contains("results")) %>%
  gather(target, predictions, contains("predictions")) %>%
  unnest(predictions) %>%
  group_by(trait, parameter, quantile, pEnv, rep, target) %>%
  summarize(accuracy = cor(value, predicted.value), rmse = sqrt(mean((value - predicted.value)^2)), .groups = "drop")

# Plot
a <- 0.10

# Summarize for plotting
genotype_stability_extreme_samples_pov_acc_summ <- genotype_stability_extreme_samples_pov_acc %>%
  group_by(trait, parameter, quantile, pEnv, target) %>%
  summarize_at(vars(accuracy), list(mean = mean, upper = ~quantile(., 1 - a/2, na.rm = T), lower = ~quantile(., a/2, na.rm = T)), 
               na.rm = TRUE) %>%
  ungroup() %>%
  mutate(group = paste0(pEnv, target, sep  = ".")) 


genotype_stability_extreme_samples_pov_acc_summ %>%
  ggplot(aes(x = pEnv, y = mean, color = target)) +
  # geom_ribbon(aes(ymin = lower, ymax = upper, fill = target), alpha = 0.2) +
  geom_line() +
  facet_grid(trait + parameter ~ quantile)

# New colors for different parameters
parameter_colors_use <- set_names(neyhart_palette("umn1")[3:5], names(coef_replace1))

g_prediction_samples <- genotype_stability_extreme_samples_pov_acc_summ %>%
  filter(target == "predictions_test", parameter != "sd_d") %>%
  mutate(parameter = factor(parameter, levels = names(coef_replace1))) %>%
  ggplot(aes(x = pEnv, y = mean, color = parameter)) +
  geom_hline(data = filter(genotype_stability_samples_pov_acc_summ, pEnv == 0.5,target == "predictions_test",
                           trait %in% genotype_stability_extreme_samples_pov_acc_summ$trait, parameter != "sd_d"),
             aes(yintercept = mean, color = parameter, linetype = "Average accuracy when using a\nrandom 50% of environments")) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = parameter), alpha = 0.2) +
  geom_line() +
  scale_x_continuous(name = "Proportion of environments sampled from the quantile", breaks = pretty) +
  scale_y_continuous(name = "Between-generation prediction accuracy", breaks = pretty) +
  scale_fill_manual(values = parameter_colors_use, guide = "none") +
  scale_color_manual(values = parameter_colors_use, labels = coef_replace1, name = NULL) +
  scale_linetype_manual(name = NULL, values = 2) +
  facet_grid(trait ~ quantile, switch = "y", labeller = labeller(trait = str_add_space, quantile = c("lower" = "lower 50%", "upper" = "upper 50%", "middle" = "middle 50%", "tails" = "upper/lower 25%"))) +
  theme_presentation(base_size = 10) +
  theme(strip.placement = "outside", legend.position = "top")

ggsave(filename = "figureS4_predictions_sample_environment_extremes.jpg", plot = g_prediction_samples,
       height = 6, width = 8, dpi = 1000, path = fig_dir)



# Save results ------------------------------------------------------------


save("genotype_stability_pov_acc", "genotype_stability_pov_results", "genotype_stability_kfold_acc", "genotype_stability_kfold_results",
     "genotype_stability_samples_pov_results", "genotype_stability_samples_pov_acc", "genotype_stability_extreme_samples_pov_results",
     "genotype_stability_extreme_samples_pov_acc",
     file = file.path(result_dir, "genomewide_prediction_results.RData"))

