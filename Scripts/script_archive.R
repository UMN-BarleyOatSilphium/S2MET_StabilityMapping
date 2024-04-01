## Archive code - not used for manuscript

## GWAS

# Plot
g_qq <- qq_gwas_geno_mean_adj %>% 
  ggplot(aes(x = neg_log10_p_exp, y = p_value_neg_log10, col = model)) + 
  geom_ribbon(aes(x = neg_log10_p_exp, ymin = ci_upper, ymax = ci_lower), fill = "grey75", inherit.aes = FALSE) +
  geom_abline(slope = 1, intercept = 0) + 
  geom_point() + 
  facet_wrap(~ trait, ncol = 2) +
  scale_color_discrete(name = "Model") +
  labs(title = "QQ Plot") +
  ylab(expression(-log[10](italic(p)))) +
  xlab(expression(Expected~-log[10](italic(p)))) +
  theme_bw() +
  theme(legend.position = c(0.8, 0.2))


save_file <- file.path(fig_dir, "gwas_pheno_genotype_mean_qq.jpg")
ggsave(filename = save_file, plot = g_qq, width = 5, height = 5)


## Quality control
## QQ plot
qq_gwas_geno_mean_adj <- gwas_pheno_mean_fw_adj %>% 
  filter(coef == "g") %>%
  arrange(trait, model, p_value) %>%
  group_by(trait, model) %>%
  mutate(p_exp = ppoints(n = n()), # Generate p values under null of no associations
         neg_log10_p_exp = -log10(p_exp), # Convert to -log10(p)
         # Add a confidence interval based on the beta distribution (assumes independence of tests)
         ci_lower = -log10(qbeta(p = (alpha / 2), shape1 = seq(n()), rev(seq(n())))),
         ci_upper = -log10(qbeta(p = 1 - (alpha / 2), shape1 = seq(n()), rev(seq(n()))))) %>%
  select(trait, model, marker, neg_log10_p_exp, p_value_neg_log10, ci_lower, ci_upper)



## Manhattan plot
g_man <- gwas_pheno_mean_fw_adj %>%
  filter(coef == "g") %>%
  mutate(chrom = as.factor(chrom)) %>% 
  # ggplot(aes(x = pos / 1000000, y = neg_log10_p_adj, group = chrom, col = chrom)) + 
  ggplot(aes(x = pos / 1000000, y = q_value_neg_log10, group = chrom, col = chrom)) +
  facet_grid(trait + model ~ chrom, switch = "x", scales = "free", space = "free_x") +
  g_add

save_file <- file.path(fig_dir, "gwas_pheno_genotype_mean_manhattan.jpg")
ggsave(filename = save_file, plot = g_man, width = 9, height = 12)




## Only consider the G models
g_man_Gmodels <- gwas_pheno_mean_fw_adj %>%
  filter(coef == "g") %>%
  filter(model %in% c("G", "QG")) %>%
  mutate(chrom = as.factor(chrom)) %>%
  # ggplot(aes(x = pos / 1000000, y = neg_log10_p_adj, group = chrom, col = chrom)) + 
  ggplot(aes(x = pos / 1000000, y = q_value_neg_log10, group = chrom, col = chrom)) + 
  facet_grid(trait + model ~ chrom, switch = "x", scales = "free", space = "free_x") +
  g_add

save_file <- file.path(fig_dir, "gwas_pheno_genotype_mean_manhattan_Gmodels.jpg")
ggsave(filename = save_file, plot = g_man_Gmodels, width = 9, height = 10)







## Quality control
## QQ plot
qq_gwas_pheno_stab_adj <- gwas_pheno_mean_fw_adj %>% 
  filter(coef != "g") %>%
  arrange(trait, coef, model, p_value) %>%
  group_by(trait, coef, model) %>%
  mutate(p_exp = ppoints(n = n()), # Generate p values under null of no associations
         neg_log10_p_exp = -log10(p_exp), # Convert to -log10(p)
         # Add a confidence interval based on the beta distribution (assumes independence of tests)
         ci_lower = -log10(qbeta(p = (alpha / 2), shape1 = seq(n()), rev(seq(n())))),
         ci_upper = -log10(qbeta(p = 1 - (alpha / 2), shape1 = seq(n()), rev(seq(n()))))) %>%
  select(trait, coef, model, marker, neg_log10_p_exp, p_value_neg_log10, ci_lower, ci_upper)

# Plot
g_qq <- qq_gwas_pheno_stab_adj %>% 
  ggplot(aes(x = neg_log10_p_exp, y = p_value_neg_log10, col = model)) + 
  geom_ribbon(aes(x = neg_log10_p_exp, ymin = ci_upper, ymax = ci_lower), fill = "grey75", inherit.aes = FALSE) +
  geom_abline(slope = 1, intercept = 0) + 
  geom_point() + 
  facet_grid(trait ~ coef) +
  scale_color_discrete(name = "Model") +
  labs(title = "QQ Plot") +
  ylab(expression(-log[10](italic(p)))) +
  xlab(expression(Expected~-log[10](italic(p)))) +
  theme_bw()

save_file <- file.path(fig_dir, "gwas_pheno_stability_qq.jpg")
ggsave(filename = save_file, plot = g_qq, width = 5, height = 5)




g_man_tp <- gwas_pheno_mean_fw_adj %>%
  filter(coef != "g") %>%
  mutate(chrom = as.factor(chrom)) %>% 
  # ggplot(aes(x = pos / 1000000, y = neg_log10_p_adj, group = chrom, col = chrom)) + 
  ggplot(aes(x = pos / 1000000, y = q_value_neg_log10, group = chrom, col = chrom)) +
  facet_grid(trait + coef + model ~ chrom, switch = "x", scales = "free", space = "free_x") +
  g_add

save_file <- file.path(fig_dir, "gwas_pheno_fw_manhattan.jpg")
ggsave(filename = save_file, plot = g_man_tp, width = 9, height = 12)


## Only consider "G" models
g_man_tp_Gmodels <- gwas_pheno_mean_fw_adj %>%
  filter(model %in% c("G", "QG"), coef != "g") %>%
  mutate(chrom = as.factor(chrom),
         coef = if_else(coef == "b", "Linear", "Non-Linear")) %>% 
  # ggplot(aes(x = pos / 1000000, y = neg_log10_p_adj, group = chrom, col = chrom)) + 
  ggplot(aes(x = pos / 1000000, y = q_value_neg_log10, group = chrom, col = chrom)) +
  facet_grid(trait + coef + model ~ chrom, switch = "x", scales = "free", space = "free_x") +
  g_add

save_file <- file.path(fig_dir, "gwas_pheno_fw_manhattan_Gmodels.jpg")
ggsave(filename = save_file, plot = g_man_tp_Gmodels, width = 9, height = 12)




## Combine the mean and stability manhattan plots into a single graph



# Model to use
model_toplot <- "G"

# Combine mean and stability into a single df.
# Convert the stability -log10(q) values into their opposites (i.e. *-1)
gwas_geno_mean_use <- gwas_pheno_mean_fw_adj %>%
  ungroup() %>%
  filter(model == "QG", coef == "g")  %>%
  select(marker:trait, contains("neg_log10_fdr"), Mean = q_value_neg_log10)

gwas_pheno_fw_use <- gwas_pheno_mean_fw_adj %>%
  ungroup() %>%
  filter(model == "QG", coef != "g") %>%
  mutate(coef = if_else(coef == "b", "linear", "non-linear")) %>%
  select(marker:trait, contains("neg_log10_fdr"), coef, Stability = q_value_neg_log10)

# Extract the palette for use in the manhattan plot
palette_use <- umn_palette("Secondary_Tier2")[c(3,8,4,9)]
odd_chrom <- seq(1, 7, by = 2)

gwas_use <-  gwas_pheno_fw_use %>%
  filter(coef == "linear") %>%
  full_join(gwas_geno_mean_use, .) %>%
  gather(test, q_value_neg_log10, -marker:-neg_log10_fdr10, -coef) %>%
  mutate_at(vars(q_value_neg_log10, contains("neg_log10_fdr")), funs(if_else(test == "Stability", . * -1, .))) %>%
  mutate(color = case_when(chrom %in% odd_chrom & test == "Mean" ~ "palette1",
                           !chrom %in% odd_chrom & test == "Mean" ~ "palette2",
                           chrom %in% odd_chrom & test == "Stability" ~ "palette3",
                           !chrom %in% odd_chrom & test == "Stability" ~ "palette4"),
         chrom = as.factor(chrom))

# color vector for ggplot
color_value <- set_names(palette_use, unique(gwas_use$color))

# Create the combo plot
## Create the common modifier
g_mod <- list(
  geom_hline(yintercept = 0),
  geom_point(aes(col = color)),
  geom_hline(aes(yintercept = neg_log10_fdr05, lty = "FDR 05%")),
  geom_hline(aes(yintercept = neg_log10_fdr10, lty = "FDR 10%")),
  facet_grid(trait ~ chrom, switch = "x", scales = "free_x", space = "free_x"),
  scale_linetype_discrete(guide = FALSE),
  scale_color_manual(values = color_value, labels = c("Mean", "Mean", "Stability", "Stability"),
                     guide = FALSE),
  # ylab("-log(q)"),
  # xlab("Position (Mbp)"),
  # labs(title = "Genomewide Assocation Analyses"),
  theme_manhattan()
)


## Plot GrainYield
g_mean_fw_gy <- gwas_use %>%
  filter(trait == "GrainYield") %>%
  ggplot(aes(x = pos / 1000000, y = q_value_neg_log10)) +
  ylim(c(-1.5, 1.5)) +
  g_mod +
  theme(axis.text.x = element_blank(),
        strip.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank())

## Plot HeadingDate
g_mean_fw_hd <- gwas_use %>%
  filter(trait == "HeadingDate") %>%
  ggplot(aes(x = pos / 1000000, y = q_value_neg_log10)) +
  ylim(c(-4, 4)) +
  g_mod +
  ylab(expression(-log[10](italic(q)))) +
  theme(axis.text.x = element_blank(),
        strip.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

## Plot PlantHeight
g_mean_fw_ph <- gwas_use %>%
  filter(trait == "PlantHeight") %>%
  ggplot(aes(x = pos / 1000000, y = q_value_neg_log10)) +
  ylim(c(-4.5, 4.5)) +
  g_mod +
  xlab("Position (Mbp)") +
  theme(axis.title.y = element_blank())

# Create grobs
grobs <- list(g_mean_fw_gy, g_mean_fw_hd, g_mean_fw_ph) %>%
  map(ggplotGrob)

save_file <- file.path(fig_dir, "gwas_pheno_mean_and_linear_stability.jpg")
ggsave(filename = save_file, plot = do.call("rbind", grobs), width = 9, height = 6)




## Plot non-linear stability
gwas_use <-  gwas_pheno_fw_use %>%
  filter(coef == "non-linear") %>%
  full_join(gwas_geno_mean_use, .) %>%
  gather(test, q_value_neg_log10, -marker:-neg_log10_fdr10, -coef) %>%
  mutate_at(vars(q_value_neg_log10, contains("neg_log10_fdr")), funs(if_else(test == "Stability", . * -1, .))) %>%
  mutate(color = case_when(chrom %in% odd_chrom & test == "Mean" ~ "palette1",
                           !chrom %in% odd_chrom & test == "Mean" ~ "palette2",
                           chrom %in% odd_chrom & test == "Stability" ~ "palette3",
                           !chrom %in% odd_chrom & test == "Stability" ~ "palette4"),
         chrom = as.factor(chrom))

# Create the combo plot
## Create the common modifier
g_mod <- list(
  geom_hline(yintercept = 0),
  geom_point(aes(col = color)),
  geom_hline(aes(yintercept = neg_log10_fdr05, lty = "FDR 05%")),
  geom_hline(aes(yintercept = neg_log10_fdr10, lty = "FDR 10%")),
  facet_grid(trait ~ chrom, switch = "x", scales = "free_x", space = "free_x"),
  scale_linetype_discrete(guide = FALSE),
  scale_color_manual(values = color_value, labels = c("Mean", "Mean", "Stability", "Stability"),
                     guide = FALSE),
  # ylab("-log(q)"),
  # xlab("Position (Mbp)"),
  # labs(title = "Genomewide Assocation Analyses"),
  theme_manhattan()
)


## Plot GrainYield
g_mean_fw_gy <- gwas_use %>%
  filter(trait == "GrainYield") %>%
  ggplot(aes(x = pos / 1000000, y = q_value_neg_log10)) +
  ylim(c(-1.5, 1.5)) +
  g_mod +
  theme(axis.text.x = element_blank(),
        strip.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_blank())

## Plot HeadingDate
g_mean_fw_hd <- gwas_use %>%
  filter(trait == "HeadingDate") %>%
  ggplot(aes(x = pos / 1000000, y = q_value_neg_log10)) +
  ylim(c(-7, 7)) +
  g_mod +
  ylab(expression(-log[10](italic(q)))) +
  theme(axis.text.x = element_blank(),
        strip.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

## Plot PlantHeight
g_mean_fw_ph <- gwas_use %>%
  filter(trait == "PlantHeight") %>%
  ggplot(aes(x = pos /  1000000, y = q_value_neg_log10)) +
  ylim(c(-3.5, 3.5)) +
  g_mod +
  xlab("Position (Mbp)") +
  theme(axis.title.y = element_blank())

# Create grobs
grobs <- list(g_mean_fw_gy, g_mean_fw_hd, g_mean_fw_ph) %>%
  map(ggplotGrob)

save_file <- file.path(fig_dir, "gwas_pheno_mean_and_nonlinear_stability.jpg")
ggsave(filename = save_file, plot = do.call("rbind", grobs), width = 9, height = 6)





### Marker Effect Stability


## Plot the distribution of linear and non-linear stability estimates for markers

# Get the distinct observations of the stability measures
S2_MET_marker_eff_pheno_fw_uniq <- S2_MET_marker_eff_pheno_fw %>% 
  group_by(trait, marker) %>% 
  mutate(mean_marker_effect = mean(mar_effect)) %>% 
  ungroup() %>%
  select(trait, marker:pos, b_std_error, df, stability_term, estimate, mean_marker_effect) %>% 
  distinct()

# Plot histograms of the regression coefficients and MSEs
g_mar_fw_dens <- S2_MET_marker_eff_pheno_fw_uniq %>%
  ggplot(aes(x = estimate)) + 
  geom_density(aes(fill = "blue")) + 
  facet_wrap(trait ~ stability_term, ncol = 2, scales = "free") +
  xlab("Estimate") +
  labs(title = "Marker Effect Stability Measures") +
  scale_fill_brewer(palette = "Set2", guide = FALSE) +
  theme_bw() +
  theme(axis.title.y = element_blank())

save_file <- file.path(fig_dir, "marker_stability_estimate_distriubtions.jpg")
ggsave(filename = save_file, plot = g_mar_fw_dens, height = 6, width = 5)


# Calculate the mean effect per marker and plot that against the regression coefficient
g_mar_eff_stab <- S2_MET_marker_eff_pheno_fw_uniq %>% 
  filter(stability_term == "b") %>% 
  ggplot(aes(x = mean_marker_effect, y = estimate)) + 
  geom_point() + 
  facet_grid(~ trait, scales = "free_x") +
  theme_bw()

save_file <- file.path(fig_dir, "marker_mean_effect_and_stability.jpg")
ggsave(filename = save_file, plot = g_mar_eff_stab, height = 5, width = 7)




## Plot a genomic map of marker effect stability

## Determine outlier markers using a 2.5% empirical threshold

# What should be the significance level
alpha <- 0.05

# For each trait, calculate empirical thresholds for significance
S2_MET_marker_eff_pheno_fw_sig <- S2_MET_marker_eff_pheno_fw_uniq %>%
  filter(stability_term == "b") %>% 
  group_by(trait) %>% 
  # mutate(estimate = scale(estimate)) %>%
  mutate(lower_perc = quantile(estimate, alpha / 2), 
         upper_perc = quantile(estimate, 1 - (alpha / 2))) %>%
  ungroup()

# Plot the linear stability over genomic distance
g_mar_eff_fw_b <- S2_MET_marker_eff_pheno_fw_sig %>% 
  ggplot(aes(x = pos / 1000000, y = estimate)) + 
  geom_abline(slope = 0, intercept = 0) +
  # Threshold lines
  geom_hline(aes(yintercept = lower_perc), lty = 2) +
  geom_hline(aes(yintercept = upper_perc), lty = 2) +
  geom_point() + 
  facet_grid(trait ~ chrom, space = "free_x", scales = "free", switch = "x") +
  scale_color_gradient2() +
  ylab("Stability Coefficient") +
  xlab("Position (Mbp)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.spacing.x = unit(x = 0, units = "in")) +
  labs(title = "Marker Linear Stability/Sensitivity")

# Save
save_file <- file.path(fig_dir, "marker_linear_stability.jpg")
ggsave(filename = save_file, plot = g_mar_eff_fw_b, height = 7, width = 9)


## Plot non-linear stability
# For each trait, calculate empirical thresholds for significance
S2_MET_marker_eff_pheno_fw_delta_sig <- S2_MET_marker_eff_pheno_fw_uniq %>%
  filter(stability_term == "delta") %>% 
  group_by(trait) %>% 
  # mutate(estimate = scale(estimate)) %>%
  mutate(upper_perc = quantile(estimate, 1 - (alpha / 2))) %>%
  ungroup()

g_mar_eff_fw_delta <- S2_MET_marker_eff_pheno_fw_delta_sig %>% 
  ggplot(aes(x = pos / 1000000, y = estimate)) + 
  geom_hline(aes(yintercept = upper_perc), lty = 2) +
  geom_point() + 
  facet_grid(trait ~ chrom, space = "free_x", scales = "free", switch = "x") +
  ylab("Stability Coefficient") +
  xlab("Position (Mbp)") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.spacing.x = unit(x = 0, units = "in"))  +
  labs(title = "Marker Non-Linear Stability/Sensitivity")

save_file <- file.path(fig_dir, "marker_nonlinear_stability.jpg")
ggsave(filename = save_file, plot = g_mar_eff_fw_delta, height = 7, width = 9)



## Use the slope estimates and standard errors to perform some hypothesis testing

## For the slope estimate, use a t-test to determine if the slope coefficient is 
## significantly greater than 0 (i.e. plastic) or signficantly less than 0 (i.e. stable)

## What happens if we test the marker slopes via a t-test?
mar_fw_linear_adj <- S2_MET_marker_eff_pheno_fw_uniq  %>%
  filter(stability_term == "b") %>%
  group_by(trait) %>%
  mutate(t_stat = estimate / b_std_error,
         p_not_zero = 2 * pt(abs(t_stat), df, lower.tail = FALSE),
         p_gt_zero = pt(t_stat, df, lower.tail = FALSE),
         p_lt_zero = 1 - p_gt_zero) %>%
  # Correct for multiple testing
  mutate_at(vars(starts_with("p_")), funs(qval = qvalue(.)$qvalue, local_fdr = qvalue(.)$lfdr)) %>%
  mutate_at(vars(contains("qval")), funs(neg_log = -log10(.)))

# Combine the stability vs sensitivity p-values
# gtlt = greater than / less than
mar_fw_linear_adj_gtlt <- mar_fw_linear_adj %>%
  # mutate(., p_lt_zero_qval_neg_log = -1 * p_lt_zero_qval_neg_log) %>%
  select(trait:pos, p_gt_zero_qval_neg_log, p_lt_zero_qval_neg_log, p_gt_zero_local_fdr, p_lt_zero_local_fdr) %>% 
  # Two sets of gather to tidy the dat
  gather(test_type, neg_log_q, p_gt_zero_qval_neg_log, p_lt_zero_qval_neg_log) %>% 
  gather(fdr_test, local_fdr, p_gt_zero_local_fdr, p_lt_zero_local_fdr) %>%
  filter((fdr_test == "p_gt_zero_local_fdr" & test_type == "p_gt_zero_qval_neg_log") |
           (fdr_test == "p_lt_zero_local_fdr" & test_type == "p_lt_zero_qval_neg_log")) %>%
  select(-fdr_test) %>%
  mutate(test_type = ifelse(test_type == "p_gt_zero_qval_neg_log", "plastic", "stable"),
         neg_log10_fdr05 = -log10(0.05),
         neg_log10_fdr10 = -log10(0.10))

# Filter out q-values for significance, then take the top q-value for each chromosome
top_qs <- mar_fw_linear_adj_gtlt %>%
  filter(neg_log_q >= neg_log10_fdr05) %>% 
  group_by(trait, test_type, chrom) %>% 
  top_n(1, neg_log_q) %>%
  mutate(local_fdr = str_c("fdr: ", round(local_fdr, 3)))

# Plots for both gt and lt
g_mar_fw_negpos_plots <- mar_fw_linear_adj_gtlt %>%
  # filter(trait == "HeadingDate") %>%
  mutate(., chrom = as.factor(chrom)) %>%
  ggplot(aes(x = pos / 1000000, y = neg_log_q)) +
  geom_point(aes(group = chrom, col = chrom)) + 
  # geom_text(data = top_qs, aes(x = pos / 1000000, y = neg_log_q, label = local_fdr), 
  #           size = 1.5, nudge_x = -150, check_overlap = T) +
  geom_hline(aes(yintercept = neg_log10_fdr05, lty = "FDR 05%")) +
  geom_hline(aes(yintercept = neg_log10_fdr10, lty = "FDR 10%")) +
  facet_grid(trait + test_type ~ chrom, switch = "x", scales = "free", space = "free_x") +
  scale_color_manual(guide = FALSE,
                     values = set_names(rep(c("black", "grey75"), length.out = 7), seq(1, 7))) +
  scale_shape_manual(values = c(16,1), labels = c("Sensitivity", "Stability"), name = "Test") +
  scale_linetype_discrete(guide = FALSE) +
  ylab("-log(q)") +
  xlab("Position (Mbp)") +
  labs(title = "Genomewide Test for Marker Linear Stability") +
  theme_bw() +
  theme(panel.spacing.x = unit(x = 0, units = "cm"),
        panel.spacing.y = unit(x = 0.5, units = "cm"),
        panel.border = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))


save_file <- file.path(fig_dir, "marker_linear_stability_gtle_hyptest.jpg")
ggsave(filename = save_file, plot = g_mar_fw_negpos_plots, height = 9, width = 9)


## Plot just the test of the slope being different from zero
mar_fw_linear_adj_nz <- mar_fw_linear_adj %>%
  mutate(., neg_log10_fdr05 = -log10(0.05), neg_log10_fdr10 = -log10(0.10))

g_mar_fw_nz_plots <- mar_fw_linear_adj_nz %>%
  mutate(., chrom = as.factor(chrom)) %>%
  ggplot(aes(x = pos / 1000000, y = p_not_zero_qval_neg_log)) +
  geom_point() +
  geom_hline(aes(yintercept = neg_log10_fdr05, lty = "FDR 05%")) +
  geom_hline(aes(yintercept = neg_log10_fdr10, lty = "FDR 10%")) +
  facet_grid(trait ~ chrom, switch = "x", scales = "free", space = "free_x") +
  scale_color_manual(guide = FALSE,
                     values = set_names(rep(c("black", "grey75"), length.out = 7), seq(1, 7))) +
  scale_linetype_discrete(guide = FALSE) +
  ylab("-log(p)") +
  xlab("Position (Mbp)") +
  # labs(title = "Genomewide Assocation Analysis for Phenotypic Stability",
  #      caption = "n = 183 lines used to calculate stability coefficients and n = 175 lines used\nfor association analysis. The bold line indicates the genomewide FDR threshold at 5%;\nthe dashed line at 10%.") +
  theme_bw() +
  theme(panel.spacing.x = unit(x = 0, units = "cm"),
        panel.spacing.y = unit(x = 0.5, units = "cm"),
        panel.border = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

save_file <- file.path(fig_dir, "marker_linear_stability_nz_hyptest.jpg")
ggsave(filename = save_file, plot = g_mar_fw_nz_plots, height = 7, width = 9)


# Save the adjusted p-value results from the marker effect stability analysis
save_file <- file.path(result_dir, "S2MET_pheno_mar_eff_fw_results.RData")
save("S2_MET_marker_eff_pheno_fw", "mar_fw_linear_adj_gtlt", file = save_file)



## Filter outliers and determine windows

## We will continue only with the TP data and for the slope of the regression curve. 
## We define stable marker outliers as those below the lower threshold and sensitive
## markers as those above the upper threshold

# Determine the number of stable or sensitive markers for each trait
mar_eff_fw_b_tp <- S2_MET_marker_eff_pheno_fw_sig$pop_tp %>% 
  select(trait, marker:cM_pos, estimate:upper_perc) %>% 
  mutate(outlier = case_when(estimate >= upper_perc ~ "sensitive", 
                             estimate <= lower_perc ~ "stable", 
                             TRUE ~ "not_outlier")) %>%
  distinct()

# Filter the outliers
mar_stability_outliers <- mar_eff_fw_b_tp %>%
  filter(outlier != "not_outlier")

# Summarize
mar_stability_outliers %>% 
  group_by(trait, chrom, outlier) %>% 
  summarize(n_outliers = n())




## For each trait, chromosome, and outlier type, create different windows from each
## outlier SNP to cluster other outlier SNPS

# Vector of window sizes (in bp)
wind_size_list <- c(50000, 100000, 500000, 1000000, 50000000)

mar_stability_outliers_window <- mar_stability_outliers %>% 
  split(list(.$trait, .$chrom, .$outlier)) %>% 
  .[map_dbl(., nrow) > 0] %>%
  list(., map(., function(df) {
    # Iterate over window sizes
    groups <- wind_size_list %>% 
      map(function(wind_size) {
        # Create an IRanges object
        df_iranges <- IRanges(start = df$pos - wind_size, end = df$pos + wind_size)
        
        # Reduce and find unique groups
        df_iranges %>% 
          reduce() %>% 
          findOverlaps(query = ., subject = df_iranges) %>% 
          queryHits() %>%
          str_c("chr", df$chrom, "_", .) })
    
    # Return a data.frame
    groups %>% 
      set_names(str_c("window_", wind_size_list)) %>% 
      as_data_frame() })) %>%
  pmap_df(~bind_cols(.x, .y))


















## Annotate GWAS results
## 
## This script will include code for doing the following:
## 1. Look at the LD pattern surrounding significant associations
## 2. Find overlap between associations for the mean and stability
## 3. Find overlap between associations detected in our analysis versus
## results from previous analyses
## 4. Generate examples of overlap and nearby gene annotations.


## 

# Load some packages
library(tidyverse)
library(readxl)
library(GenomicRanges)
library(neyhart)
library(cowplot)

# Directory containing the project repository
repo_dir <- getwd()

# Project and other directories
source(file.path(repo_dir, "source.R"))

# Load the genome annotation data
load(file.path(result_dir, "snp_genome_annotation.RData"))
# Read in the significant markers identified in our analysis
load(file.path(result_dir, "gwas_adjusted_significant_results.RData"))

# Add a window of 5 Mbp to either side of the significant markers
window <- 5e6

# Create a marker matrix M
M <- S2TP_imputed_multi_genos_mat

# Data frame of SNP information (pos, chrom, etc)
snp_info <- S2TP_imputed_multi_genos_hmp %>% 
  select(marker = rs, chrom, pos, cM_pos)


### First examine the LD structure around significant markers

# Gather significant markers and create a GRange object
gwas_sig_mar_grange <- gwas_mlmm_marker_info %>%
  # Add twice the window to complement the other overlapping procedures
  mutate(start = pos - (2 * window), end = pos + (2 * window)) %>%
  makeGRangesFromDataFrame(., keep.extra.columns = TRUE)

# Gather all markers and create a Grange object
all_marker_grange <- snp_info %>% 
  select(marker, chrom, start = pos, cM_pos) %>% 
  mutate(end = start) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T)

# For each significant marker, find the markers within the interval
gwas_sig_mar_overlap <- findOverlaps(query = gwas_sig_mar_grange, subject = all_marker_grange) %>%
  as.data.frame() %>% 
  split(.$queryHit) %>% 
  map(~list(query_marker = unique(gwas_sig_mar_grange[.$queryHits]$marker), 
            subject_marker = all_marker_grange[.$subjectHits]$marker))

# Now for each query marker, find the LD between that marker and the subject markers
gwas_sig_mar_overlap_LD <- gwas_sig_mar_overlap %>% 
  map(function(mars) M[,c(mars$query_marker, mars$subject_marker)] %>% 
  {cor(.)^2} %>% .[mars$query_marker, mars$subject_marker,drop = FALSE]) %>%
  map_df(~data.frame(query_marker = row.names(.), subject_marker = colnames(.), LD = c(.)))

# Get the positions of these SNPS
gwas_sig_mar_overlap_LD_pos <- gwas_sig_mar_overlap_LD %>% 
  left_join(., snp_info, by = c("query_marker" = "marker")) %>% 
  select(chrom, query_marker, query_pos = pos, subject_marker, LD) %>% 
  left_join(., snp_info, by = c("subject_marker" = "marker")) %>% 
  select(chrom = chrom.x, query_marker:subject_marker, subject_marker_pos = pos, LD)

# Filter out the same markers
gwas_sig_mar_overlap_LD_pos1 <- gwas_sig_mar_overlap_LD_pos %>%
  filter(query_marker != subject_marker)

# Plot
g_sig_mar_LD <- gwas_sig_mar_overlap_LD_pos1 %>%
  filter(query_marker == unique(.$query_marker)) %>% 
  ggplot(aes(x = subject_marker_pos / 1e6, y = LD)) +
  geom_point() + 
  geom_vline(aes(xintercept = query_pos / 1e6)) + 
  facet_wrap(~ query_marker, scales = "free_x")



### Next determine the significant markers for the mean per se that overlap with
### those for stability

# Split the significant marker Grange by trait
# If the trait does not have at least one coefficient, return NULL
# Reset the window so it isn't twice the size
gwas_sig_mar_grange_split <- gwas_sig_mar_grange %>% 
  `start<-`(., value = .$pos - window) %>% 
  `end<-`(., value = .$pos + window) %>%
  as.data.frame() %>%
  split(.$trait)

gwas_sig_mar_grange_split1 <- gwas_sig_mar_grange_split[map_dbl(gwas_sig_mar_grange_split, ~n_distinct(.$coef)) > 1]

## Now we are only looking at heading date and plant height
# Split by coefficient
gwas_sig_mar_grange_split2 <- gwas_sig_mar_grange_split1 %>%
  map(~split(., .$coef) %>% map(makeGRangesFromDataFrame, keep.extra.columns = T))

# For each trait, iterate over the non-mean per se list and find overlaps with the mean per se list
gwas_sig_mar_grange_overlap <- gwas_sig_mar_grange_split2 %>%
  map_df(function(trait_gr) {
    query_coef <- "g"
    subject_coef <- setdiff(names(trait_gr), query_coef)
    
    # Iterate over the stability coefficients
    hits <- map(subject_coef, ~mergeByOverlaps(query = trait_gr[[query_coef]], subject = trait_gr[[.]]))
    
    # Use the hits to find the overlapping markers
    hits %>% 
      map_df(~as.list(.) %>% as.data.frame(.) %>% 
               select(trait, chrom = trait_gr..query_coef...seqnames, mps_marker = marker,
                      mps_marker_pos = pos, mps_alpha = alpha, stab_coef = coef.1,
                      stab_marker = marker.1,stab_marker_pos = pos.1, stab_alpha = alpha.1))
  })

## Which is the minor allele?
af <- colMeans(M + 1) / 2

## What is the LD between these marker pairs?
gwas_sig_mar_grange_overlap %>% 
  left_join(x = ., y = select(gwas_sig_mar_overlap_LD_pos1, query_marker, subject_marker, LD), 
            by = c("mps_marker" = "query_marker", "stab_marker" = "subject_marker")) %>%
  mutate(mps_marker_af = af[mps_marker], stab_marker_af = af[stab_marker])







### Check the significant markers for overlap with loci previously associated with
### the mean per se of agronomic traits

# Read in the QTL metadata files that were downloaded from T3
qtl_meta <- map_df(list.files(data_dir, pattern = "meta", full.names = TRUE), read_csv)
# Read in association data from Pauli2014 and Wang2012
bcap_association <- read_csv(file.path(data_dir, "BCAP_association_qtl.csv")) %>%
  mutate(position = parse_number(position)) %>%
  select(trait, marker, chrom = chromosome, position, gene:feature, reference) %>%
  filter(trait %in% c("GrainYield", "HeadingDate", "PlantHeight"))


# Rename the traits
qtl_meta_df <- qtl_meta %>%
  mutate(trait = case_when(trait == "grain yield" ~ "GrainYield",
                           trait == "heading date" ~ "HeadingDate",
                           trait == "plant height" ~ "PlantHeight"))

# Remove unmapped QTL and change the chromosome name
qtl_meta_use <- qtl_meta_df %>% 
  filter(!chromosome %in% c("chrUNK", "chrUn")) %>% 
  mutate(chrom = as.integer(parse_number(chromosome))) %>%
  select(trait, marker, chrom, position, gene, feature) %>%
  arrange(trait, chrom, position)

# Combine data
qtl_meta_use1 <- bind_rows(qtl_meta_use, bcap_association)



## Create GRanges for each marker data.frame
qtl_meta_grange <- qtl_meta_use1 %>%
  split(.$trait) %>%
  map(~makeGRangesFromDataFrame(., keep.extra.columns = TRUE, start.field = "position", 
                                end.field = "position"))



# Use the window from above to find previously-identified markers that coincide with
# those discovered in our anaylsis.

# Create a GRange per trait
gwas_mlmm_grange <- gwas_mlmm_marker_info %>% 
  mutate(start = pos - window, end = pos + window) %>%
  split(.$trait) %>%
  map(~makeGRangesFromDataFrame(., keep.extra.columns = TRUE))

## For each trait, see how many of the markers detected in our analysis overlap
## with those previously detected in the germplasm?
gwas_markers_overlap <- list(gwas_mlmm_grange, qtl_meta_grange) %>%
  pmap_df(function(gwas, t3) {
    overlaps <- findOverlaps(query = gwas, subject = t3)
    
    # Add the overlapping markers to the gwas marker grange
    overlap_data <- bind_cols(
      as.data.frame(gwas[queryHits(overlaps)]),
      as.data.frame(t3[subjectHits(overlaps)])
    )
    
    # Bind with the markers that did not overlaps
    non_overlap_data <- as.data.frame(gwas[setdiff(seq_along(gwas), queryHits(overlaps))])
    
    # Sort and remove superfluous columns
    bind_rows(overlap_data, non_overlap_data) %>%
      select(trait, coef, marker, chrom = seqnames, pos, start, end, alpha, qvalue, R_sqr_snp, af, 
             qtl_marker = marker1, qtl_pos = start1, reference) %>% 
      arrange(trait, coef, chrom, pos)
    
    
  })

## Convert af to MAF and parse the coefficients
gwas_markers_overlap_toprint <- gwas_markers_overlap %>%
  mutate(coef = str_replace_all(coef, coef_replace),
         MAF = pmin(af, 1 - af)) %>% 
  select(trait:qvalue, MAF, names(.)) %>%
  arrange(trait, coef, chrom, pos)


## Write to a table
save_file <- file.path(fig_dir, "gwas_significant_associations.csv")
write_csv(x = gwas_markers_overlap_toprint, path = save_file)

## Re-annotate by hand
# Read in the annotated excel file
gwas_markers_overlap <- file.path(fig_dir, "gwas_significant_associations_annotated.xlsx") %>%
  read_excel(na = c("NA")) %>% 
  rename_all(str_to_lower) %>% 
  dplyr::rename(coef = character, chrom = chromosome)



## Summarize this table into a condensed version for the main manuscript
gwas_markers_overlap_summ <- gwas_markers_overlap %>%
  group_by(trait, coef, marker) %>% 
  summarize(overlapped = any(!is.na(coincidentcharacter)), ref = list(reference)) %>%
  summarize(n_overlap = sum(overlapped), n_mar = n(), ref = list(ref)) %>% 
  ungroup() %>%
  mutate(ref = map_chr(ref, ~unlist(.) %>% str_split(";") %>% unlist() %>% 
                         str_trim() %>% unique() %>% na.omit() %>% str_c(collapse = "; "))) %>%
  dplyr::rename(Trait = trait, Character = coef, `Coincident Markers` = n_overlap,
                `Total Markers` = n_mar, References = ref)

## Write this
save_file <- file.path(fig_dir, "gwas_sig_marker_annotated_summary.csv")
write_csv(x = gwas_markers_overlap_summ, path = save_file)




### Summary statistics
# Number of associations per trait-coef pair
# Also show the distinct chromosomes
gwas_markers_overlap %>%
  group_by(trait, coef) %>%
  nest(marker, chrom) %>%
  mutate(n_association = map_int(data, ~n_distinct(.$marker)),
         chrom = map(data, ~unique(.$chrom))) %>%
  select(-data) %>% 
  as.data.frame()

## How many markers overlapped with previous associations?
## Calculate by trait and coefficient
gwas_markers_overlap %>%
  group_by(trait, coef, marker) %>% 
  summarize(overlapped = any(!is.na(coincidentcharacter))) %>% 
  summarize(n_overlap = sum(overlapped), prop_overlapped = n_overlap / n())

## Quick notes
# 3 / 8 HD linear stability QTL overlapped with genotype mean QTL
# 5 / 8 HD genotype mean QTL overlapped with stability QTL

## What is the mean amount of variation explained by markers influencing stability
## that overlapped with mean per se versus not overlapped?
gwas_markers_overlap_rsqr <- gwas_markers_overlap %>% 
  group_by(trait, coef, marker) %>%
  mutate(overlapped = any(!is.na(coincidentcharacter))) %>% 
  distinct(trait, coef, marker, alpha, `r^2`, overlapped) %>% 
  group_by(trait, coef, overlapped) %>% 
  summarize(n_mar = n(), R_sqr_snp_mean = mean(`r^2`),
            R_sqr_snp = list(`r^2`))


gwas_markers_overlap_rsqr_test <-gwas_markers_overlap_rsqr %>% 
  do(rank_sum_test = {
    df <- .
    if (nrow(df) > 1) {
      wilcox.test(x = df$R_sqr_snp[[1]], y = df$R_sqr_snp[[2]])
    } else {
      NA
    }
  })





### Visualize the overlaps
# Using the overlaps in the above analysis, visualize the overlaps by comparing plots

# Function to label y axis
y_lab <- function(x) format(x, nsmall = 1)

## First show the manhattan plots with annotation
# Manhanttan plot modifier
g_mod_man <- list(
  geom_point(),
  geom_hline(yintercept = -log10(0.05), lty = 2),
  # geom_hline(aes(yintercept = neg_log10_fdr10, lty = "FDR 10%")),
  scale_color_manual(values = color, guide = FALSE),
  scale_y_continuous(labels = y_lab),
  ylab(expression(-log[10](italic(q)))),
  xlab("Position (Mbp)"),
  theme_bw(),
  theme_manhattan(),
  theme(panel.border = element_blank())
)

## Replace the data in the 'gwas_pheno_mean_fw_adj' data frame with the data in the 
## 'gwas_mlmm_marker_info' data frame if markers are in the 'gwas_mlmm_marker_info'
## data.frame
gwas_marker_sig <- gwas_mlmm_marker_info %>%
  select(trait, coef, marker, chrom) %>%
  mutate(color = "Bl")

gwas_marker_notsig <- gwas_pheno_mean_fw_adj %>%
  select(trait, coef, marker, chrom) %>%
  dplyr::setdiff(., select(gwas_marker_sig, -color)) %>%
  mutate(color = if_else(chrom %in% seq(1, 7, 2), "B", "G"))

## Replace data
gwas_results_toplot <- bind_rows(
  left_join(gwas_marker_notsig, gwas_pheno_mean_fw_adj) %>% 
    select(trait, coef, marker, chrom, pos, beta, pvalue, qvalue, color),
  left_join(gwas_marker_sig, gwas_mlmm_marker_info) %>% 
    select(trait, coef, marker, chrom, pos, beta = alpha, pvalue, qvalue, color)
) %>%
  mutate(qvalue_neg_log10 = -log10(qvalue),
         coef = str_replace_all(coef, coef_replace))



## Iterate over traits to create the completed manhattan plot
for (tr in unique(gwas_results_toplot$trait)) {
  
  # Subset the qtl meta data
  qtl_meta_use1_tr <- qtl_meta_use1 %>% 
    filter(trait == tr) %>% 
    select(trait:chrom, pos = position) %>% 
    distinct()
  
  ## First construct the manhattan plot
  g_man <- gwas_results_toplot %>%
    filter(trait == tr) %>%
    ggplot(aes(x = pos / 1000000, y = qvalue_neg_log10, group = chrom, col = color)) + 
    facet_grid(coef ~ chrom, switch = "x", scales = "free_x", space = "free_x") +
    g_mod_man + 
    ylab(expression(-log[10](italic(q))~'for single-trait')) +
    labs(title = tr) +
    theme(strip.background.x = element_blank(),
          strip.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank())
  
  # Next construct the pleiotropy plot
  g_man_plei <- gwas_pheno_mean_fw_plei_toplot %>%
    filter(trait == tr) %>%
    ggplot(aes(x = pos / 1000000, y = -log10(qvalue), color = color)) +
    facet_grid(stab_coef ~ chrom, switch = "x", scales = "free_x", space = "free_x") +
    g_mod_man +
    ylab(expression(atop(-log[10](italic(q))~'for pleiotropy', 'with Genotype Mean'))) +
    theme(strip.background.x = element_blank(),
          strip.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank())
  
  # Next construct the annotation plot
  g_ann <- snp_info %>%
    ggplot(aes(x = pos / 1e6, y = 0)) +  # Need to use the snp data to get the x-axes aligned
    geom_point(color = "white") + 
    # geom_text_repel(size = 2, ) +
    geom_segment(data = qtl_meta_use1_tr, lwd = 5, 
                 mapping = aes(x = (pos - 1e6) / 1e6, xend = (pos + 1e6) / 1e6, y = 0, yend = 0)) +
    xlab("Position (Mbp)") + 
    facet_grid(. ~ chrom, switch = "x", scales = "free", space = "free_x") +
    theme_manhattan() +
    theme(panel.border = element_blank(),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank())
  
  ## Combine
  g_man_ann <- plot_grid(g_man, g_man_plei, g_ann, ncol = 1, rel_heights = c(1, 2/3, 0.25), 
                         align = "hv", axis = "lr", labels = c("A", "B", "C"))
  
  # Save
  save_file <- file.path(fig_dir, str_c("gwas_manhattan_complete_annotated_", tr, ".jpg"))
  ggsave(filename = save_file, plot = g_man_ann, width = 9, height = 10, dpi = 1000)
  
}





### To show more closely the coincidence of heading date mean per se associations
### and linear stability associations, create a combine plot with the results for
### mean per se, the results for linear stability, and the pleiotropy results


# Function to label y axis
y_lab <- function(x) format(x, nsmall = 1)

# Create a common plot modifier
g_mod_man <- list(
  geom_point(),
  scale_y_continuous(labels = y_lab),
  scale_color_manual(values = color, guide = FALSE),
  ylab(expression(-log[10](italic(q)))),
  theme_bw(),
  theme(panel.border = element_blank(),
        panel.grid = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y = element_line())
)


### Example for heading date
tr <- "HeadingDate"
chr <- 5
start <- 580e6
end <- 620e6

# Extract the results for mean per se and linear stability
gwas_hd_example <- gwas_pheno_mean_fw_adj %>% 
  filter(trait == tr, coef %in% c("g", "b"), 
         chrom == chr, between(pos, start, end))

# Extract the gene annotations
gene_annotations <- snp_info_grange_overlap$genes %>% 
  as.data.frame() %>% 
  select(chrom = 1, gene_start = ..start, gene_end = ..end, gene_id) %>% 
  mutate(chrom = parse_number(chrom)) %>% 
  filter(chrom == chr, between(gene_start, start, end), 
         between(gene_end, start, end)) %>% 
  distinct()

## Data.frame for signficant BOPA hits for heading date
qtl_meta_use_hd <- qtl_meta_use1 %>% 
  filter(trait == tr, chrom == chr, between(position, start, end)) %>%
  mutate(gene_start = position, gene_end = position) %>%
  select(chrom, gene_start, gene_end, gene_id = marker)



# First plot the results for mean per se
g_hd_example_mean <- gwas_hd_example %>% 
  filter(coef == "g") %>% 
  ggplot(aes(x = pos / 1000000, y = -log10(qvalue))) + 
  g_mod_man +
  ylab(expression(-log[10](italic(q))~'for mean'~italic(per~se))) +
  labs(title = tr) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), units = "cm"),
        plot.title = element_text(hjust = 0.5))

# Second plot the results for linear stability
g_hd_example_fw <- gwas_hd_example %>% 
  filter(coef == "b") %>% ggplot(aes(x = pos / 1000000, y = log10(qvalue))) + 
  g_mod_man +
  ylab(expression(log[10](italic(q))~'for linear stability'))  +
  theme(plot.margin = unit(c(0, 0, 0, 0), units = "cm"))


g_hd_annotations <- gene_annotations %>%
  # g_hd_annotations <- gene_annotations_use %>% 
  ggplot(aes(x = gene_start / 1000000, xend = gene_end / 1000000, y = 0, yend = 0)) + 
  geom_segment(lwd = 5, color = "white") + # Change color to white to hide
  # Add VRN-H1
  geom_segment(data = vrn1_data,
               aes(x = (gene_start - 100000) / 1000000, xend = (gene_end + 100000) / 1000000,
                   y = 1, yend = 1), inherit.aes = T, lwd = 5, color = color[2]) +
  geom_text(data = vrn1_data, aes(label = gene_id, y = 1), inherit.aes = T, vjust = 2, fontface = "italic") +
  # Add significant QTL
  geom_segment(data = qtl_meta_use_hd,
               aes(x = (gene_start - 100000) / 1000000, xend = (gene_end + 100000) / 1000000,
                   y = 0.5, yend = 0.5), inherit.aes = T, lwd = 5, color = color[1]) +
  # geom_text(data = qtl_meta_use_hd, aes(label = gene_id, y = 0.5), inherit.aes = T, vjust = 2, 
  #           check_overlap = T, size = 3, hjust = "inward") +
  xlab("Position on Chromosome 5 (Mbp)") +
  ylim(c(0, 1)) +
  theme_bw() +
  theme(axis.line.x = element_line(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(c(0, 0.5, 0, 0.5), units = "cm"))

# 
## Subset the pleiotropy results
g_hd_example_plei <- gwas_pheno_mean_fw_plei_toplot %>%
  filter(trait == tr, stab_coef == "Linear Stability", chrom == chr,
         between(pos, start, end)) %>%
  ggplot(aes(x = pos / 1000000, y = -log10(qvalue))) +
  g_mod_man +
  ylab(expression(-log[10](italic(q))~'for pleiotropy')) +
  theme(panel.border = element_rect(fill = "transparent"),
        axis.text.x = element_text(),
        axis.ticks.x = element_line(),
        axis.line.x = element_line())


# Combine
hd_example_combine <- plot_grid(
  g_hd_example_mean, g_hd_annotations, g_hd_example_fw, 
  g_hd_example_plei,
  ncol = 1, align = "hv", 
  rel_heights = c(1, 0.4, 1, 0.4))

## Save
save_file <- file.path(fig_dir, "hd_gwas_annotation_example.jpg")
ggsave(filename = save_file, plot = hd_example_combine, width = 6, height = 10, dpi = 1000)





## Plot the example for height
tr <- "PlantHeight"
chr <- 6
start <- 0
end <- 20e6

# Extract the results for mean per se and linear stability
gwas_ph_example <- gwas_pheno_mean_fw_adj %>% 
  filter(trait == tr, coef %in% c("g", "b"), chrom == chr, between(pos, start, end))

## Data.frame for signficant BOPA hits for height
qtl_meta_use_ph <- qtl_meta_use1 %>% 
  filter(trait == tr, chrom == chr, between(position, start, end)) %>%
  mutate(gene_start = position, gene_end = position) %>%
  select(chrom, gene_start, gene_end, gene_id = marker)


# First plot the results for mean per se
g_ph_example_mean <- gwas_ph_example %>% 
  filter(coef == "g") %>% 
  ggplot(aes(x = pos / 1000000, y = -log10(qvalue))) + 
  g_mod_man +
  ylab(expression(-log[10](italic(q))~'for mean'~italic(per~se))) +
  labs(title = tr) + 
  theme(plot.margin = unit(c(0, 0, 0, 0), units = "cm"),
        plot.title = element_text(hjust = 0.5))

# Second plot the results for linear stability
g_ph_example_fw <- gwas_ph_example %>% 
  filter(coef == "b") %>% ggplot(aes(x = pos / 1000000, y = log10(qvalue))) + 
  g_mod_man +
  ylab(expression(log[10](italic(q))~'for linear stability'))  +
  theme(plot.margin = unit(c(0, 0, 0, 0), units = "cm"))

gene_annotations <- snp_info_grange_overlap$genes %>% 
  as.data.frame() %>% 
  select(chrom = 1, gene_start = ..start, gene_end = ..end, gene_id) %>% 
  mutate(chrom = parse_number(chrom)) %>% 
  filter(chrom == chr, between(gene_start, start, end), 
         between(gene_end, start, end)) %>% 
  distinct()

# Add annotation
g_ph_annotations <- gene_annotations %>% 
  ggplot(aes(x = gene_start / 1000000, xend = gene_end / 1000000, y = 0, yend = 0)) + 
  geom_segment(lwd = 5, color = "white") + # Change color to white to hide
  # Add significant QTL
  geom_segment(data = qtl_meta_use_ph,
               aes(x = (gene_start - 100000) / 1000000, xend = (gene_end + 100000) / 1000000,
                   y = 0.5, yend = 0.5), inherit.aes = T, lwd = 5, color = color[1]) +
  # geom_text(data = qtl_meta_use_ph, aes(label = gene_id, y = 0.5), inherit.aes = T, vjust = 2, 
  #           check_overlap = T, size = 3, hjust = "inward") +
  xlab("Position on Chromosome 7 (Mbp)") +
  ylim(c(0, 1)) +
  theme_bw() +
  theme(axis.line.x = element_line(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        plot.margin = unit(c(0, 0.5, 0, 0.5), units = "cm"))

## Subset the pleiotropy results
g_ph_example_plei <- gwas_pheno_mean_fw_plei_toplot %>%
  filter(trait == tr, stab_coef == "Linear Stability", chrom == chr,
         between(pos, start, end)) %>%
  ggplot(aes(x = pos / 1000000, y = -log10(qvalue))) +
  g_mod_man +
  ylab(expression(-log[10](italic(q))~'for pleiotropy')) +
  theme(panel.border = element_rect(fill = "transparent"),
        axis.text.x = element_text(),
        axis.ticks.x = element_line(),
        axis.line.x = element_line())


# Combine
ph_example_combine <- plot_grid(
  g_ph_example_mean, g_ph_annotations, g_ph_example_fw, 
  g_ph_example_plei,
  ncol = 1, align = "hv", 
  rel_heights = c(1, 0.4, 1, 0.4))

## Save
save_file <- file.path(fig_dir, "ph_gwas_annotation_example.jpg")
ggsave(filename = save_file, plot = ph_example_combine, width = 6, height = 10, dpi = 1000)



## Combine the heading date and plant height plots
hd_ph_example_combine <- plot_grid(hd_example_combine, ph_example_combine, ncol = 2,
                                   align = "hv", labels = c("A", "B"))

# Save
save_file <- file.path(fig_dir, "hd_ph_gwas_annotation_example.jpg")
ggsave(filename = save_file, plot = hd_ph_example_combine, width = 10, height = 10, dpi = 1000)




