## S2MET Mapping
## 
## Stability prediction analysis
## 
## This notebook will outline genome-wide predictions of stability/sensitivity using the S2MET data

library(tidyverse)
library(readxl)
library(lme4)
library(broom)
library(stringr)

repo_dir <- getwd()

# Project and other directories
source(file.path(repo_dir, "source.R"))

# Load the results
load(file.path(result_dir, "S2MET_stability_crossv_results.RData"))


## Cross-Validation of Training Set using all markers or marker subsets
## The first results use relationship matrices created using all markers in each subgroup

# Edit the results for plotting
cv_results_toplot <- cv_results %>%
  mutate(coef = str_replace_all(coef, coef_replace),
         coef = str_replace_all(coef, pattern = " ", replacement = "\n"),
         coef = parse_factor(coef, levels = c("Genotype\nMean", "Linear\nStability", "Non-Linear\nStability")))

## Only plot the results of using all markers
g_cv_acc_all <- cv_results_toplot %>%
  filter(significance == "all_markers") %>%
  ggplot(aes(x = coef, y = acc)) +
  geom_boxplot(position = "dodge") +
  ylab("Prediction Accuracy") +
  xlab("") +
  ylim(c(-0.4, 1)) +
  # labs(caption = "Results of 100 iterations of 7-fold cross-validation") +
  facet_grid(~trait) +
  # theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_bw()

# Save this
save_file <- file.path(fig_dir, "cv_results_all.jpg")
ggsave(filename = save_file, plot = g_cv_acc_all, width = 8, height = 5)

# Calculate average accuracy
(cv_results_avg_acc <- cv_results %>% 
  filter(significance == "all_markers") %>%
  group_by(trait, coef) %>% 
  summarize(mean_acc = mean(acc)) %>%
  spread(coef, mean_acc))





## Results of using marker samples of the same size

# First, take the average accuracy of the 10 replications of each sample
cv_results_samples_avg1 <- cv_results_samples_tidy %>% 
  group_by(trait, coef, rep, significance) %>% 
  summarize(acc = mean(acc.rep_acc)) %>%
  ungroup()

# Create a df for plotting
cv_results_samples_toplot <- cv_results_samples_avg1 %>%
  mutate(coef = str_replace_all(coef, coef_replace),
         coef = str_replace_all(coef, pattern = " ", replacement = "\n"),
         coef = parse_factor(coef, levels = c("Genotype\nMean", "Linear\nStability", "Non-Linear\nStability")),
         significance = str_replace_all(significance, c("all" = "All", "average" = "Average", "sensitive" = "Plastic",
                                                        "stable" = "Stable")))

# Plot boxplots
g_cv_acc_samples <- cv_results_samples_toplot %>% 
  ggplot(aes(x = coef, y = acc, color = significance)) + 
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


(cv_results_samples_avg2 <- cv_results_samples_avg1 %>% 
  group_by(trait, coef, significance) %>% 
  summarize(mean_acc = mean(acc)) %>% 
  spread(significance, mean_acc))

# Calculate the percent change over random markers
(cv_results_samples_perc <- cv_results_samples_avg2 %>% 
  mutate_at(vars(average:stable), funs((. - all) / all)))




