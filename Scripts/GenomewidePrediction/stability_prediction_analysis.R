## S2MET Mapping
## 
## Stability prediction analysis
## 
## This notebook will outline genome-wide predictions of stability/sensitivity using the S2MET data

repo_dir <- getwd()
# Project and other directories
source(file.path(repo_dir, "source.R"))

# Load the results
load(file.path(result_dir, "stability_crossv_results.RData"))
# Load the stability results
load(file.path(result_dir, "pheno_mean_fw_results.RData"))
# Load the vp prediction results
load(file.path(result_dir, "vp_stability_prediction.RData"))


# Reassign the marker matrix
M <- S2TP_imputed_multi_genos_mat






## Cross-Validation of Training Set using all markers or marker subsets
## The first results use relationship matrices created using all markers in each subgroup

## Summarize the accuracy for each rep
cv_results <- cv_results_samples_tidy %>% 
  group_by(trait, coef, rep, significance) %>% 
  summarize(acc = mean(acc.rep_acc)) %>%
  ungroup()



# Edit the results for plotting
cv_results_toplot <- cv_results %>%
  mutate(coef = str_replace_all(coef, coef_replace),
         coef = str_replace_all(coef, pattern = " ", replacement = "\n"),
         coef = parse_factor(coef, levels = c("Genotype\nMean", "Linear\nStability", "Non-Linear\nStability")))

## Only plot the results of using all markers
g_cv_acc_all <- cv_results_toplot %>%
  filter(significance == "all") %>%
  ggplot(aes(x = coef, y = acc)) +
  geom_boxplot(position = "dodge") +
  ylab("Prediction Accuracy") +
  xlab("") +
  ylim(c(-0.4, 1)) +
  labs(title = "Cross-Validation Using All Markers",
       caption = str_c("Results of ", n_distinct(cv_results_toplot$rep), " iterations of 7-fold cross-validation")) +
  facet_grid(~trait) +
  theme_bw()

# Save this
ggsave(filename = "cv_results_all_markers.jpg", plot = g_cv_acc_all, path = fig_dir,
       width = 6, height = 4, dpi = 1000)




## Results of using marker samples of the same size
# Significance level
alpha <- 0.05

# Colors for the different marker groups
colors <- c("black", umn_palette(n = 4)[4]) %>%
  set_names(c("Stable", "Plastic"))


## Calculate the mean, sd, and 95% confidence interval
cv_results_samples_summ <- cv_results %>% 
  group_by(trait, coef, significance) %>% 
  summarize_at(vars(acc), funs(mean_acc = mean(.), sd_acc = sd(.), n = n())) %>% 
  mutate(se_acc = sd_acc / sqrt(n), conf = qt(p = 1 - (alpha / 2), df = n - 1) * se_acc, 
         lower = mean_acc - conf, 
         upper = mean_acc + conf) %>%
  ungroup()


# Create a df for plotting
# cv_results_samples_toplot <- cv_results_samples_summ %>%
cv_results_samples_toplot <- cv_results %>%
  mutate(coef = str_replace_all(coef, coef_replace),
         coef = str_replace_all(coef, pattern = " ", replacement = "\n"),
         coef = parse_factor(coef, levels = c("Genotype\nMean", "Linear\nStability", "Non-Linear\nStability")),
         significance = str_replace_all(significance, c("all" = "All", "plastic" = "Plastic", "stable" = "Stable", 
                                                        "mixed" = "Mixed")))

# Plot boxplots
g_cv_acc_samples <- cv_results_samples_toplot %>% 
  ggplot(aes(x = coef, y = acc, fill = significance)) + 
  geom_boxplot(position = "dodge") + 
  ylab("Prediction Accuracy") +
  xlab("") +
  ylim(c(-0.4, 1)) + 
  # labs(caption = "Results of 100 iterations of 7-fold cross-validation using constant marker size (232).") +
  scale_color_discrete(name = "Marker Group") +
  facet_grid(~trait) + 
  theme_bw()


# Save this
save_file <- file.path(fig_dir, "cv_results_samples.jpg")
ggsave(filename = save_file, plot = g_cv_acc_samples, width = 8, height = 5)


### Plot for the poster

g_cv_acc_samples_poster <- cv_results_samples_toplot %>% 
  mutate(coef = as.character(coef),
         coef = ifelse(coef == "Linear\nStability", "Phenotypic\nPlasticity", coef)) %>%
  filter(coef != "Non-Linear\nStability") %>%
  ggplot(aes(x = coef, y = mean_acc, fill = significance, group = significance)) + 
  # geom_point(position = position_dodge(0.9)) + 
  geom_col(position = "dodge") +
  # geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.9), width = 0.5) +
  geom_errorbar(aes(ymin = mean_acc - sd_acc, ymax = mean_acc + sd_acc),
                position = position_dodge(0.9), width = 0.5) +
  ylab("Prediction Accuracy") +
  xlab("") +
  # ylim(c(-0.4, 1)) + 
  ylim(c(0, 1)) +
  scale_fill_manual(values = c("All" = "grey75", colors), name = "Marker Type") +
  facet_wrap(~trait, ncol = 2) + 
  labs(caption = "n = 232 plastic and stable markers\n100 samples of n = 232 markers drawn from all markers and from the average markers.") +
  theme_poster() +
  theme(legend.position = c(0.75, 0.25),
        title = element_text(size = 10))

ggsave(filename = "cv_results_samples_poster.jpg", plot = g_cv_acc_samples_poster, 
       height = 9, width = 8, path = fig_dir, dpi = 1000)




(cv_results_samples_avg2 <- cv_results_samples_avg1 %>% 
  group_by(trait, coef, significance) %>% 
  summarize(mean_acc = mean(acc)) %>% 
  spread(significance, mean_acc))

# Calculate the percent change over random markers
(cv_results_samples_perc <- cv_results_samples_avg2 %>% 
  mutate_at(vars(average:stable), funs((. - all) / all)))



## VP prediction results
# Summarize all marker VP prediction
all_marker_prediction_results$vp_prediction_all_data %>%
  summarize(acc = cor(value, pred_value)) %>%
  spread(coef, acc)
    
# trait            b     g log_delta
# 1 GrainYield  0.288  0.553    0.0532
# 2 HeadingDate 0.473  0.535    0.198 
# 3 PlantHeight 0.0914 0.614    0.332 


# Summarize the LOYO results
vp_prediction_loyo_summ <- all_marker_prediction_results$vp_prediction_loyo %>% 
  bind_rows() %>% 
  group_by(year, trait, coef) %>%
  summarize(acc = cor(value, pred_value, use = "complete.obs"))

# Plot
vp_prediction_loyo_summ %>% 
  ungroup() %>% 
  mutate(year = as.factor(year)) %>% 
  ggplot(aes(x = coef, y = acc, fill = year, group = year)) + 
  geom_col(position = "dodge") + 
  facet_wrap(~trait, ncol = 2)

# Plot
g_loyo_acc <- vp_prediction_loyo_summ %>% 
  ungroup() %>% 
  mutate(year = as.factor(year)) %>% 
  ggplot(aes(x = coef, y = acc, fill = year, group = year)) + 
  geom_col(position = "dodge") + 
  scale_fill_discrete(name = NULL) +
  facet_wrap(~trait, ncol = 3) +
  theme_bw() + 
  theme(legend.position = c(0.06, 0.85), legend.text = element_text(size = 6))

ggsave(filename = "vp_prediction_loyo.jpg", plot = g_loyo_acc, path = fig_dir, 
       height = 4, width = 8, dpi = 1000)


## Summarize the random environment results
vp_prediction_rand_env_summ <- all_marker_prediction_results$vp_prediction_rand_env %>% 
  list(., seq_along(.)) %>%
  pmap_df(~mutate(.x, iter = .y)) %>%
  group_by(iter, trait, coef) %>%
  summarize(acc = cor(value, pred_value, use = "complete.obs"))

# Plot
g_rand_env_acc <- vp_prediction_rand_env_summ %>% 
  ungroup() %>% 
  ggplot(aes(x = coef, y = acc, fill = coef)) + 
  geom_boxplot() + 
  scale_fill_discrete(name = NULL) +
  facet_wrap(~trait, ncol = 3) +
  theme_bw() + 
  theme(legend.position = c(0.06, 0.85), legend.text = element_text(size = 6))

ggsave(filename = "vp_prediction_rand_env.jpg", plot = g_rand_env_acc, path = fig_dir, 
       height = 4, width = 8, dpi = 1000)






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












