## S2MET Mapping
## 
## Stability prediction analysis
## 
## This notebook will outline genome-wide predictions of stability/sensitivity using the S2MET data

repo_dir <- getwd()
# Project and other directories
source(file.path(repo_dir, "source.R"))

# Other packages
library(cowplot)

# Load the results
load(file.path(result_dir, "genomewide_prediction_results.RData"))
# Load the stability results
load(file.path(result_dir, "pheno_mean_fw_results.RData"))



## Color scheme
colors_use <- set_names(umn_palette(2)[3:5], coef_replace)




## Cross-Validation of Training Set using all markers or marker subsets
## The first results use relationship matrices created using all markers in each subgroup

# Compare correlations
cv_results_alldata %>% 
  summarize_at(vars(acc), funs(mean = mean(.), sd = sd(.)))

# trait       coef         mean     sd
# 1 GrainYield  b         0.431   0.0689
# 2 GrainYield  g         0.448   0.0621
# 3 GrainYield  log_delta 0.138   0.0882
# 4 HeadingDate b         0.427   0.0825
# 5 HeadingDate g         0.858   0.0226
# 6 HeadingDate log_delta 0.0487  0.0993
# 7 PlantHeight b         0.403   0.0853
# 8 PlantHeight g         0.517   0.0688
# 9 PlantHeight log_delta 0.00602 0.0964

# Boxplot
g_cv_results_box <- cv_results_alldata %>% 
  ungroup() %>%
  mutate(coef = str_replace_all(coef, coef_replace)) %>%
  ggplot(aes(x = trait, y = acc, fill = coef)) + 
  geom_boxplot(position = "dodge") + 
  scale_fill_manual(values = colors_use, name = NULL) +
  ylab("Prediction accuracy") +
  ylim(c(-0.3, 1)) +
  theme_pnas() + 
  theme(axis.title.x = element_blank(), legend.position = c(0.80, 0.90))

# Save
ggsave(filename = "cv_all_boxplot.jpg", plot = g_cv_results_box, path = fig_dir,
       height = 9, width = 8, units = "cm", dpi = 1000)





## Cross-validation for g and b when sampling environments
cv_results_test_env <- cv_results_env %>%
  select(trait:pred_gb) %>%
  unnest() %>%
  rename_at(vars(contains("acc")), funs(str_replace(., "acc_pred_", ""))) %>%
  gather(coef, acc, g:log_delta)

# Compare correlations
cv_results_test_env %>% 
  group_by(trait, coef) %>%
  summarize_at(vars(acc), funs(mean = mean(.), sd = sd(.)))

# trait       coef   mean     sd
# 1 GrainYield  b     0.179 0.181 
# 2 GrainYield  g     0.386 0.107 
# 3 HeadingDate b     0.126 0.154 
# 4 HeadingDate g     0.836 0.0364
# 5 PlantHeight b     0.269 0.136 
# 6 PlantHeight g     0.465 0.0822


# Boxplot
g_cv_results_box <- cv_results_test_env %>% 
  ungroup() %>%
  mutate(coef = str_replace_all(coef, coef_replace)) %>%
  ggplot(aes(x = trait, y = acc, fill = coef)) + 
  geom_boxplot(position = "dodge") + 
  scale_fill_manual(values = colors_use, name = NULL) +
  ylab("Prediction accuracy") +
  ylim(c(-0.5, 1)) +
  theme_pnas() + 
  theme(axis.title.x = element_blank(), legend.position = c(0.80, 0.90))

# Save
ggsave(filename = "cv_test_env_boxplot.jpg", plot = g_cv_results_box, path = fig_dir,
       height = 9, width = 8, units = "cm", dpi = 1000)


## Combine data
cv_results_coef <- bind_rows(
  mutate(cv_results_alldata, test = "AllEnvironments"),
  mutate(cv_results_test_env, test = "EnvironmentSampling")
)


# Boxplot
g_cv_results_box <- cv_results_coef %>% 
  ungroup() %>%
  mutate(coef = str_replace_all(coef, coef_replace)) %>%
  ggplot(aes(x = trait, y = acc, fill = coef)) + 
  geom_boxplot(position = "dodge") + 
  scale_fill_manual(values = colors_use, name = NULL) +
  ylab("Prediction accuracy") +
  ylim(c(-0.5, 1)) +
  facet_grid(~test) +
  theme_pnas() + 
  theme(axis.title.x = element_blank(), legend.position = c(0.13, 0.13), legend.margin = margin(),
        panel.spacing.x = unit(1, "lines"))


# Save
ggsave(filename = "cv_coef_boxplot.jpg", plot = g_cv_results_box, path = fig_dir,
       height = 9, width = 11.4, units = "cm", dpi = 1000)





## Look at environment-specific prediction accuracy
# Format the heritability results
env_herit_use <- stage_one_data %>%
  ungroup() %>% 
  filter(!str_detect(trial, "S2C1F4")) %>% 
  distinct(trait, environment, heritability) %>%
  # Filter out low heritabilities
  filter(heritability >= 0.15)


cv_results_env_spec <- cv_results_env %>%
  select(-pred_gb) %>%
  unnest()

## Add heritability and standardize the predictions
cv_results_env_spec1 <- cv_results_env_spec %>% 
  ungroup() %>%
  right_join(., env_herit_use) %>% 
  mutate_at(vars(contains("acc")), funs(. / sqrt(heritability))) %>%
  filter(!is.na(acc_gb))

# Color scheme for prediction types
colors_use_pred <- set_names(umn_palette(4)[2:3], c(acc_gb = "Mean + Stability", acc_g = "Mean"))
pred_type_replace <- set_names(names(colors_use_pred), c("acc_gb", "acc_g"))


# Plot the accuracy for each environment
g_cv_env_obs <- cv_results_env_spec1 %>%
  gather(pred_type, acc, acc_gb, acc_g) %>%
  mutate(pred_type = str_replace_all(pred_type, pred_type_replace)) %>%
  ggplot(aes(x = environment, y = acc, fill = pred_type)) + 
  geom_boxplot() +
  facet_grid(trait ~ ., switch = "y", scales = "free_y") +
  scale_fill_manual(values = colors_use_pred, name = NULL) +
  ylab("Prediction accuracy") +
  theme_pnas() +
  theme(axis.text.x = element_text(angle = 90), axis.title.x = element_blank(),
        panel.spacing.y = unit(1, "lines"), legend.position = "bottom")

## Save
ggsave(filename = "cv_env_spec_obs_boxplot.jpg", plot = g_cv_env_obs, path = fig_dir,
       height = 12, width = 17.4, units = "cm", dpi = 1000)


## Simply calculate environment means over iterations
cv_results_env_summ <- cv_results_env_spec1 %>%
  group_by(trait, environment) %>%
  summarize_at(vars(acc_gb, acc_g), mean)

# Overall mean
cv_results_env_summ %>%
  summarize_at(vars(acc_gb, acc_g), mean)

# trait       acc_gb acc_g
# 1 GrainYield   0.317 0.291
# 2 HeadingDate  0.817 0.817
# 3 PlantHeight  0.394 0.386


## Compare the distributions of these means using a t test
# The test is paired because it is the same TP that is used for prediction in each iteration
cv_results_env_summ %>% 
  do(t_test = t.test(x = .$acc_gb, y = .$acc_g, paired = TRUE, conf.int = TRUE)) %>%
  ungroup() %>%
  mutate(t_test_pvalue = map_dbl(t_test, "p.value"))

# trait       t_test      t_test_pvalue
# 1 GrainYield  <S3: htest>         0.180
# 2 HeadingDate <S3: htest>         0.792
# 3 PlantHeight <S3: htest>         0.315


## Plot the distribution of mean standardized accuracy
g_cv_env_mean <- cv_results_env_summ %>%
  gather(pred_type, acc, acc_gb, acc_g) %>%
  mutate(pred_type = str_replace_all(pred_type, pred_type_replace)) %>%
  ggplot(aes(x = trait, y = acc, fill = pred_type)) + 
  geom_boxplot() +
  scale_fill_manual(values = colors_use_pred, name = NULL) +
  ylab("Standardized prediction accuracy") +
  ylim(c(-0.6, 1.2)) + 
  theme_pnas() +
  theme(axis.title.x = element_blank(), legend.position = c(0.82, 0.13))

## Save
ggsave(filename = "cv_env_mean_boxplot.jpg", plot = g_cv_env_mean, path = fig_dir,
       height = 8, width = 8.7, units = "cm", dpi = 1000)


## Combine the cross-validation plots
g_cv_combined <- plot_grid(
  g_cv_results_box + ylim(c(-0.6, 1.2)) + theme(legend.direction = "horizontal", legend.position = c(0.5, 0.95)),
  g_cv_env_mean + ylim(c(-0.6, 1.2)) + theme(legend.position = c(0.80, 0.95)),
  nrow = 1, align = "hv", axis = "tblr", labels = LETTERS[1:2], label_size = 10, rel_widths = c(1, 0.66)
)

ggsave(filename = "cv_combined.jpg", plot = g_cv_combined, path = fig_dir,
       height = 8, width = 17.8, units = "cm", dpi = 1000)











### Predictions of the VP

## Predict g and b using all the data
vp_prediction_summ <- vp_prediction_all_data %>% 
  summarize(acc = cor(value, pred_value)) %>%
  mutate(annotation = str_c("r[MG]~'= ", formatC(acc, digits = 3, format = "f"), "'"),
         coef = str_replace_all(coef, coef_replace))

vp_prediction_summ %>%
  select(-annotation) %>%
  spread(coef, acc)

# trait       `Genotype Mean` `Linear Stability` `Non-Linear Stability`
# 1 GrainYield            0.553             0.289                  0.0494
# 2 HeadingDate           0.535             0.472                  0.215 
# 3 PlantHeight           0.614             0.0944                 0.297 



## Create a plot per coef, then combine
# Prepare a df
data_toplot <- vp_prediction_all_data %>%
  ungroup() %>%
  mutate(coef = str_replace_all(coef, coef_replace),
         coef = as.factor(coef))

# Mean
g_all_data_predict_g <- data_toplot %>% 
  filter(coef == "Genotype Mean") %>%
  ggplot(aes(x = value, pred_value, color = coef, shape = coef)) + 
  geom_point(size = 0.5) + 
  geom_smooth(method = "lm", se = FALSE, lwd = 0.5) + 
  geom_text(data = filter(vp_prediction_summ, coef == "Genotype Mean"), parse = TRUE,
            aes(x = Inf, y = -Inf, label = annotation), inherit.aes = FALSE, size = 1.5, hjust = 1, vjust = -1) +
  facet_wrap(~ trait, scales = "free", strip.position = "left", ncol = 1) + 
  scale_color_manual(values = colors_use, guide = FALSE) + 
  scale_shape_discrete(guide = FALSE) + 
  ylab("Predicted value") +
  xlab("Observed value") +
  theme_pnas()

# Linear stability
g_all_data_predict_b <- data_toplot %>% 
  filter(coef == "Linear Stability") %>%
  ggplot(aes(x = value, pred_value, color = coef, shape = coef)) + 
  geom_point(size = 0.5) + 
  geom_smooth(method = "lm", se = FALSE, lwd = 0.5) + 
  geom_text(data = filter(vp_prediction_summ, coef == "Linear Stability"), parse = TRUE,
            aes(x = Inf, y = -Inf, label = annotation), inherit.aes = FALSE, size = 1.5, hjust = 1, vjust = -1) +
  facet_wrap(~ trait, scales = "free", strip.position = "left", ncol = 1) + 
  scale_color_manual(values = colors_use, name = NULL, drop = FALSE) + 
  scale_shape_discrete(guide = FALSE) + 
  ylab("Predicted value") +
  xlab("Observed value") +
  theme_pnas() +
  theme(legend.position = "bottom")

# Non-linear stability
g_all_data_predict_d <- data_toplot %>% 
  filter(coef == "Non-Linear Stability") %>%
  ggplot(aes(x = value, pred_value, color = coef, shape = coef)) + 
  geom_point(size = 0.5) + 
  geom_smooth(method = "lm", se = FALSE, lwd = 0.5) + 
  geom_text(data = filter(vp_prediction_summ, coef == "Non-Linear Stability"), parse = TRUE,
            aes(x = Inf, y = -Inf, label = annotation), inherit.aes = FALSE, size = 1.5, hjust = 1, vjust = -1) +
  facet_wrap(~ trait, scales = "free", strip.position = "left", ncol = 1) + 
  scale_color_manual(values = colors_use, guide = FALSE) + 
  scale_shape_discrete(guide = FALSE) + 
  ylab("Predicted value") +
  xlab("Observed value") +
  theme_pnas()


## Combine plots
g_vp_all_data_predict <- plot_grid(
  g_all_data_predict_g + theme(axis.title.x = element_blank(), strip.placement = "outside"),
  g_all_data_predict_b + theme(axis.title.y = element_blank(), strip.background = element_blank(),
                               strip.text = element_blank(), legend.position = "none"),
  g_all_data_predict_d + theme(axis.title = element_blank(), strip.background = element_blank(),
                               strip.text = element_blank()),
  ncol = 3, align = "hv", axis = "tb", rel_widths = c(1, 0.8, 0.8)
)

# Add legend
g_vp_all_data_predict1 <- plot_grid(
  g_vp_all_data_predict,
  get_legend(g_all_data_predict_b),
  ncol = 1, rel_heights = c(1, 0.15)
)

# Save
ggsave(filename = "vp_prediction_all_data_merge.jpg", plot = g_vp_all_data_predict1, path = fig_dir,
       height = 7, width = 8.7, units = "cm", dpi = 1000)




## Parent-offspring validation for g and b when sampling environments
vp_results_test_env <- vp_prediction_env %>%
  select(trait:pred_gb) %>%
  unnest() %>%
  rename_at(vars(contains("acc")), funs(str_replace(., "acc_pred_", ""))) %>%
  gather(coef, acc, g:b)

# Compare correlations
vp_results_test_env %>% 
  group_by(trait, coef) %>%
  summarize_at(vars(acc), funs(mean = mean(.), sd = sd(.)))

# trait       coef   mean     sd
# 1 GrainYield  b     0.0738  0.0844
# 2 GrainYield  g     0.360   0.0813
# 3 HeadingDate b     0.219   0.136 
# 4 HeadingDate g     0.453   0.0479
# 5 PlantHeight b     0.00989 0.131 
# 6 PlantHeight g     0.583   0.0275


# Boxplot
g_vp_results_box <- vp_results_test_env %>% 
  ungroup() %>%
  mutate(coef = str_replace_all(coef, coef_replace)) %>%
  ggplot(aes(x = trait, y = acc, fill = coef)) + 
  geom_boxplot(position = "dodge") + 
  scale_fill_manual(values = colors_use, name = NULL) +
  ylab("Prediction accuracy") +
  ylim(c(-0.5, 1)) +
  theme_pnas() + 
  theme(axis.title.x = element_blank(), legend.position = c(0.80, 0.90))

# Save
ggsave(filename = "cv_test_env_boxplot.jpg", plot = g_cv_results_box, path = fig_dir,
       height = 9, width = 8, units = "cm", dpi = 1000)






## Look at environment-specific prediction accuracy
# Format the heritability results
env_herit_use <- stage_one_data %>%
  ungroup() %>% 
  filter(!str_detect(trial, "S2C1F4")) %>% 
  distinct(trait, environment, heritability) %>%
  # Filter out low heritabilities
  filter(heritability >= 0.15)


vp_results_env_spec <- vp_prediction_env %>%
  select(-pred_gb) %>%
  unnest()

## Add heritability and standardize the predictions
vp_results_env_spec1 <- vp_results_env_spec %>% 
  ungroup() %>%
  right_join(., env_herit_use) %>% 
  mutate_at(vars(contains("acc")), funs(. / sqrt(heritability))) %>%
  filter(!is.na(acc_gb))

# Plot the accuracy for each environment
g_vp_env_obs <- vp_results_env_spec1 %>%
  gather(pred_type, acc, acc_gb, acc_g) %>%
  mutate(pred_type = str_replace_all(pred_type, pred_type_replace)) %>%
  ggplot(aes(x = environment, y = acc, fill = pred_type)) + 
  geom_boxplot() +
  facet_grid(trait ~ ., switch = "y", scales = "free_y") +
  scale_fill_manual(values = colors_use_pred, name = NULL) +
  ylab("Prediction accuracy") +
  theme_pnas() +
  theme(axis.text.x = element_text(angle = 90), axis.title.x = element_blank(),
        panel.spacing.y = unit(1, "lines"), legend.position = "bottom")

## Save
ggsave(filename = "vp_env_spec_obs_boxplot.jpg", plot = g_vp_env_obs, path = fig_dir,
       height = 12, width = 17.4, units = "cm", dpi = 1000)


## Simply calculate environment means over iterations
vp_results_env_summ <- vp_results_env_spec1 %>%
  group_by(trait, environment) %>%
  summarize_at(vars(acc_gb, acc_g), mean)

# Overall mean
vp_results_env_summ %>%
  summarize_at(vars(acc_gb, acc_g), mean)

# trait       acc_gb acc_g
# 1 GrainYield   0.292 0.286
# 2 HeadingDate  0.482 0.479
# 3 PlantHeight  0.555 0.560


## Compare the distributions of these means using a t test
# The test is paired because it is the same TP that is used for prediction in each iteration
vp_results_env_summ %>% 
  do(t_test = t.test(x = .$acc_gb, y = .$acc_g, paired = TRUE, conf.int = TRUE)) %>%
  ungroup() %>%
  mutate(t_test_pvalue = map_dbl(t_test, "p.value"))

# trait       t_test      t_test_pvalue
# 1 GrainYield  <S3: htest>         0.525
# 2 HeadingDate <S3: htest>         0.497
# 3 PlantHeight <S3: htest>         0.323


## Plot the distribution of mean standardized accuracy
g_vp_env_mean <- vp_results_env_summ %>%
  gather(pred_type, acc, acc_gb, acc_g) %>%
  mutate(pred_type = str_replace_all(pred_type, pred_type_replace)) %>%
  ggplot(aes(x = trait, y = acc, fill = pred_type)) + 
  geom_boxplot() +
  scale_fill_manual(values = colors_use_pred, name = NULL) +
  ylab("Standardized prediction accuracy") +
  ylim(c(-0.6, 1.2)) + 
  theme_pnas() +
  theme(axis.title.x = element_blank(), legend.position = c(0.82, 0.13))

## Save
ggsave(filename = "vp_env_mean_boxplot.jpg", plot = g_vp_env_mean, path = fig_dir,
       height = 8, width = 8.7, units = "cm", dpi = 1000)


## Combine the cross-validation plots
g_vp_combined <- plot_grid(
  g_vp_all_data_predict,
  g_vp_env_mean + ylim(c(-0.6, 1.2)) + theme(legend.position = "none"),
  nrow = 1, align = "hv", axis = "tblr", labels = LETTERS[1:2], label_size = 10, rel_widths = c(1, 0.66)
)

ggsave(filename = "vp_combined.jpg", plot = g_vp_combined, path = fig_dir,
       height = 8, width = 17.8, units = "cm", dpi = 1000)


## Combine plots
g_pred_combined <- plot_grid(
  g_cv_results_box + ylim(c(-1, 1)) + theme(legend.direction = "horizontal", legend.position = c(0.5, 0.12)),
  g_cv_env_mean + ylim(c(-0.6, 1)) + theme(legend.position = "none"),
  g_vp_all_data_predict,
  g_vp_env_mean + ylim(c(-0.6, 1)) + theme(legend.position = c(0.5, 0.99), legend.margin = margin(),
                                             legend.direction = "horizontal"),
  nrow = 2, align = "hv", axis = "tblr", labels = LETTERS[1:4], label_size = 10, rel_widths = c(1, 0.66),
  rel_heights = c(0.8, 1)
)

ggsave(filename = "prediction_combined.jpg", plot = g_pred_combined, path = fig_dir,
       height = 12, width = 17.8, units = "cm", dpi = 1000)
















## Marker subset prediction accuracy

## Random markers
rand_marker_prediction <- marker_subset_prediction_results$rand_marker_predctions %>% 
  unnest() %>% 
  group_by(trait, coef, nmar, iter) %>% 
  summarize(acc = cor(value, pred_value)) %>%
  ungroup() %>%
  mutate(nmar = parse_number(nmar)) %>% 
  group_by(trait, coef, nmar) %>%
  summarize_at(vars(acc), funs(mean = mean(.), lower = quantile(., alpha / 2), upper = quantile(., 1 - (alpha / 2)))) %>%
  ungroup()

## Plot
g_rand <- rand_marker_prediction %>% 
  ggplot(aes(x = nmar, y = mean, color = coef)) + 
  geom_point(pch = 0) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1) + 
  facet_wrap(~trait, ncol = 2) +
  theme_bw()


## Compare different marker subsets for accuracy
# Top markers
top_marker_prediction <- marker_subset_prediction_results$top_rank_markers_predictions %>% 
  list(., names(.)) %>% 
  pmap_df(~mutate(.x, nmar = .y)) %>% 
  group_by(trait, coef, nmar) %>%
  summarize(acc = cor(value, pred_value)) %>%
  ungroup() %>%
  mutate(nmar = parse_number(nmar))

## Plot
g_top <- top_marker_prediction %>% 
  ggplot(aes(x = nmar, y = acc, color = coef)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~trait, ncol = 3) +
  theme_bw()

# Add random
g_top1 <- g_top + 
  geom_point(data = rand_marker_prediction, aes(y = mean), pch = 0) +
  geom_line(data = rand_marker_prediction, aes(y = mean)) +
  geom_ribbon(data = rand_marker_prediction, aes(ymin = lower, ymax = upper, y = mean), alpha = 0.1) +
  scale_color_discrete(name = NULL) +
  theme(legend.position = c(0.90, 0.15), legend.text = element_text(size = 6))

# Save
ggsave(filename = "vp_prediction_top_marker.jpg", plot = g_top1, path = fig_dir,
       height = 4, width = 8, dpi = 1000)


# Top markers, evenly spaced
tesm_marker_prediction <- marker_subset_prediction_results$tesm_markers_predictions %>% 
  list(., names(.)) %>% 
  pmap_df(~mutate(.x, nmar = .y)) %>% 
  group_by(trait, coef, nmar) %>%
  summarize(acc = cor(value, pred_value)) %>%
  ungroup() %>%
  mutate(nmar = parse_number(nmar))

## Plot
g_tesm <- tesm_marker_prediction %>% 
  ggplot(aes(x = nmar, y = acc, color = coef)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~trait, ncol = 3) +
  theme_bw()

# Add random
g_tesm1 <- g_tesm + 
  geom_point(data = rand_marker_prediction, aes(y = mean), pch = 0) +
  geom_line(data = rand_marker_prediction, aes(y = mean)) +
  geom_ribbon(data = rand_marker_prediction, aes(ymin = lower, ymax = upper, y = mean), alpha = 0.1) +
  scale_color_discrete(name = NULL) +
  theme(legend.position = c(0.90, 0.15), legend.text = element_text(size = 6))

# Save
ggsave(filename = "vp_prediction_top_marker_evenly_spaced.jpg", plot = g_tesm1, path = fig_dir,
       height = 4, width = 8, dpi = 1000)


## Evenly-spaced markers
esm_marker_prediction <- marker_subset_prediction_results$esm_markers_predictions %>% 
  list(., names(.)) %>% 
  pmap_df(~mutate(.x, nmar = .y)) %>% 
  group_by(trait, coef, nmar) %>%
  summarize(acc = cor(value, pred_value)) %>%
  ungroup() %>%
  mutate(nmar = parse_number(nmar))

## Plot
g_esm <- esm_marker_prediction %>% 
  ggplot(aes(x = nmar, y = acc, color = coef)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~trait, ncol = 3) +
  theme_bw()

# Add random
g_esm1 <- g_esm + 
  geom_point(data = rand_marker_prediction, aes(y = mean), pch = 0) +
  geom_line(data = rand_marker_prediction, aes(y = mean)) +
  geom_ribbon(data = rand_marker_prediction, aes(ymin = lower, ymax = upper, y = mean), alpha = 0.1) +
  scale_color_discrete(name = NULL) +
  theme(legend.position = c(0.90, 0.15), legend.text = element_text(size = 6))

# Save
ggsave(filename = "vp_prediction_evenly_spaced_marker.jpg", plot = g_esm1, path = fig_dir,
       height = 4, width = 8, dpi = 1000)












