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



