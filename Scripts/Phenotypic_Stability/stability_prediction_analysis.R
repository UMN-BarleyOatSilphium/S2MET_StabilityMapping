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
# Load the GWAS results
load(file.path(result_dir, "gwas_adjusted_significant_results.RData"))

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




