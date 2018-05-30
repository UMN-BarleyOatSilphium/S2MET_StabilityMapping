## S2MET Mapping
## 
## Stability prediction analysis
## 
## This notebook will outline genome-wide predictions of stability/sensitivity using the S2MET data

library(tidyverse)
library(readxl)
library(rrBLUP)
library(lme4)
library(broom)
library(neyhart)

repo_dir <- getwd()

# Project and other directories
source(file.path(repo_dir, "source.R"))

# Load the results
load(file.path(result_dir, "S2MET_stability_crossv_results.RData"))
# Load the stability results
load(file.path(result_dir, "S2MET_pheno_mean_fw_results.RData"))
# Load the GWAS results
load(file.path(result_dir, "gwas_adjusted_significant_results.RData"))

# Reassign the marker matrix
M <- S2TP_imputed_multi_genos_mat



## Examine the marker effects between the mean per se and the stability
pheno_mean_fw_tomodel <- S2_MET_pheno_mean_fw %>%
  distinct(trait, line_name, g, b, delta) %>%
  mutate(log_delta = log(delta)) %>%
  select(-delta) %>%
  gather(coef, value, g:log_delta)

# Iterate over traits and coefficients to calculate marker effects
pheno_mean_fw_marker_eff <- pheno_mean_fw_tomodel %>%
  spread(coef, value) %>%
  group_by(trait) %>%
  do({
    df <- .
    
    # Group by coefficient and calculate marker effects
    single_coef_marker_eff <- df %>% 
      select(g, b, log_delta) %>% 
      map(~mixed.solve(y = ., Z = M)$u) %>%
      as.data.frame()
        
    ## Now fit multi-variate models
    Y <- select(df, g, b, log_delta) %>% as.matrix()
    # List of matrices
    Y_list <- list(
      g_b = Y[,1:2],
      g_delta = Y[,c(1,3)]
    )
    
    X <- model.matrix(~ 1, df)
    K <- diag(ncol(M))
    
    # Fit pairs
    fit_list <- map(Y_list, ~emmremlMultivariate(Y = t(.), X = t(X), Z = t(M), K = K))
    
    # Get the marker effect predictions
    multi_coef_marker_eff <- fit_list %>% 
      map(function(fit) fit$Gpred %>% t() %>% `colnames<-`(., row.names(fit$Vg)) %>% 
            as.data.frame(.) %>% mutate(marker = colnames(M)) %>% select(marker, names(.)))
    
    # Export a data_frame
    data_frame(single_coef = list(single_coef_marker_eff), multi_coef = list(multi_coef_marker_eff))
    
  })


# Examine the correlation in marker effects
pheno_mean_fw_marker_eff_cor <- pheno_mean_fw_marker_eff %>% 
  unnest(single_coef) %>% 
  group_by(trait) %>% 
  do(as.data.frame(cor(select(., -trait))) %>% rownames_to_column("coef"))

# Edit data to plot
pheno_mean_fw_marker_eff_toplot <- pheno_mean_fw_marker_eff %>% 
  mutate(single_coef = map(single_coef, ~rownames_to_column(., "marker"))) %>%
  unnest(single_coef) %>% 
  ungroup() %>%
  gather(coef, value, g:log_delta) %>%
  # Designate the significant markers
  full_join(., distinct(gwas_mlmm_marker_info, trait, marker, qvalue)) %>%
  mutate(significant = if_else(is.na(qvalue), "Not Significant", "Significant")) %>%
  select(trait, marker, significant, coef, value) %>%
  spread(coef, value) %>%
  rename_at(vars(b:log_delta), str_replace_all, pattern = coef_replace) %>%
  gather(coef, value, `Linear Stability`, `Non-Linear Stability`)

label_coef <- function(x) str_c(x, "\nGenotype Mean")

## Plot the marker effects
g_single_coef_plot <- pheno_mean_fw_marker_eff_toplot %>% 
  ggplot(aes(x = `Genotype Mean`, y = value, color = significant, size = significant)) + 
  geom_point() + 
  scale_color_manual(values = c('Significant' = umn_palette(name = 2)[3], "Not Significant" = "black")) +
  facet_grid(coef ~ trait, scales = "free", labeller = labeller(trait = label_coef)) +
  theme_bw() + 
  theme(axis.title = element_blank(),
        legend.position = "bottom",
        legend.title = element_blank())

# Save
save_file <- file.path(fig_dir, "single_trait_marker_effect_relationship.jpg")
ggsave(filename = save_file, plot = g_single_coef_plot, height = 6, width = 8, dpi = 1000)









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
# Significance level
alpha <- 0.05

# Colors for the different marker groups
colors <- c("black", umn_palette(n = 4)[-c(1:2)]) %>%
  set_names(c("Average", "Plastic", "Stable"))


# First, take the average accuracy of the 10 replications of each sample
cv_results_samples_avg1 <- cv_results_samples_tidy %>% 
  group_by(trait, coef, rep, significance) %>% 
  summarize(acc = mean(acc.rep_acc)) %>%
  ungroup()

## Calculate the mean, sd, and 95% confidence interval
cv_results_samples_summ <- cv_results_samples_avg1 %>% 
  group_by(trait, coef, significance) %>% 
  summarize_at(vars(acc), funs(mean_acc = mean(.), sd_acc = sd(.), n = n())) %>% 
  mutate(se_acc = sd_acc / sqrt(n), conf = qt(p = 1 - (alpha / 2), df = n - 1) * se_acc, 
         lower = mean_acc - conf, 
         upper = mean_acc + conf) %>%
  ungroup()


# Create a df for plotting
# cv_results_samples_toplot <- cv_results_samples_avg1 %>%
cv_results_samples_toplot <- cv_results_samples_summ %>%
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




