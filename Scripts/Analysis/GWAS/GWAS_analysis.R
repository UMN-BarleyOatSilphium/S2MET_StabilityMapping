## Analysis of GWAS results

# This notebook will outline the analysis of GWAS mapping procedures using the S2MET data to
# identify marker-trait associations for the mean per se of traits and phenotypic stability
# This script also includes code for generating figures related to this analysis

library(tidyverse)
library(readxl)
library(rrBLUP)
library(broom)
library(qvalue)
library(cowplot)
library(neyhart)

# Directory containing the project repository
repo_dir <- getwd()

# Project and other directories
source(file.path(repo_dir, "source.R"))

# Rename the marker matrix
M <- S2TP_imputed_multi_genos_mat

# Get the marker information
snp_info <- S2TP_imputed_multi_genos_hmp %>%
  select(marker = rs, chrom, pos, cM_pos) %>%
  # Correct the cM position for BOPA snps
  mutate(cM_pos = if_else(str_detect(marker, "^S"), cM_pos, cM_pos / 1000))

# Subset the entry list for the tp
tp_entry_list <- entry_list %>%
  filter(Line %in% tp_geno) %>%
  select(line_name = Line, program = Program)

# Significance threshold for GWAS
alpha <- 0.05


# Color scheme for manhattan plot
color <- c(set_names(umn_palette(n = 4)[3:4], "Bl", "Or"), "B" = "black", "G" = "grey75")



### Population Structure

# First examine population structure within the germplasm

# Calculate the kinship matrix
K <- A.mat(X = M, min.MAF = 0, max.missing = 1)

## Visualize pop structure
K_prcomp <- prcomp(K)

# Tidy, calculate lambda, then add program information
K_prcomp_df <- tidy(K_prcomp) %>%
  mutate(PC = str_c("PC", PC)) %>%
  dplyr::rename(line_name = row) %>%
  left_join(entry_list, by = c("line_name" = "Line"))

# Extract lambda and calculate variance explained
var_exp <- K_prcomp$sdev %>% 
  set_names(str_c("PC", seq_along(.))) %>%
  {. / sum(.)}

## Combine data.frame to plot
df_PC1_PC2 <- K_prcomp_df %>% 
  filter(PC %in% c("PC1", "PC2")) %>% 
  spread(PC, value) %>%
  mutate(x = "PC1", y = "PC2", var_expx = var_exp["PC1"], var_expy = var_exp["PC2"]) %>%
  dplyr::rename(xvalue = PC1, yvalue = PC2)

df_PC1_PC3 <- K_prcomp_df %>% 
  filter(PC %in% c("PC1", "PC3")) %>% 
  spread(PC, value) %>%
  mutate(x = "PC1", y = "PC3", var_expx = var_exp["PC1"], var_expy = var_exp["PC3"]) %>%
  dplyr::rename(xvalue = PC1, yvalue = PC3)

df_PC2_PC3 <- K_prcomp_df %>% 
  filter(PC %in% c("PC2", "PC3")) %>% 
  spread(PC, value) %>%
  mutate(x = "PC2", y = "PC3", var_expx = var_exp["PC2"], var_expy = var_exp["PC3"]) %>%
  dplyr::rename(xvalue = PC2, yvalue = PC3)

## Combine the data frames and combine the variance explained with
## the PC
df_combined <- bind_rows(df_PC1_PC2, df_PC1_PC3, df_PC2_PC3) %>%
  mutate(x1 = str_c(x, " (", 100 * round(var_expx, 3), "%)"),
         y1 = str_c(y, " (", 100 * round(var_expy, 3), "%)"))


## Plot
g_pop_str <- df_combined %>% 
  ggplot(aes(x = xvalue, y = yvalue, col = Program)) + 
  geom_point() + 
  facet_grid(y1 ~ x1, switch = "both") +
  scale_color_discrete(guide = guide_legend(title = "Breeding\nProgram")) +
  theme_bw() + 
  theme(axis.title = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank(),
        legend.position = c(0.8, 0.8))

# Save this
save_file <- file.path(fig_dir, "population_structure.jpg")
# Publication name
save_file <- file.path(fig_dir, "SupplementaryFigure1.jpg")
ggsave(filename = save_file, plot = g_pop_str, width = 5, height = 5)





## How is population structure correlated with the traits?

# Load the genotype means and FW regression results
load(file.path(result_dir, "S2MET_pheno_mean_fw_results.RData"))

# Transform the delta statistic to log-delta
S2_MET_pheno_mean_fw_trans <- S2_MET_pheno_mean_fw %>% 
  distinct(trait, line_name, g, b, delta) %>% 
  mutate(log_delta = log(delta)) %>%
  select(-delta)


# Combine the trait and PCs, duplicating the BLUEs to make summarizing smoother
trait_pop_str <- K_prcomp_df %>% 
  filter(PC %in% c("PC1", "PC2", "PC3")) %>% 
  select(line_name, program = Program, PC, eigenvalue = value) %>%
  left_join(., S2_MET_pheno_mean_fw_trans) %>%
  gather(measure, value, g:log_delta)


## Correlate the eigenvalues with the trait mean and stability, and calculate regression coefficients
## Use a bootstrap to estimate confidence interval
trait_pop_str_corr <- trait_pop_str %>% 
  group_by(trait, PC, measure) %>% 
  do(gws::boot_cor(x = .$value, y = .$eigenvalue, boot.reps = 1000)) %>%
  # Which ones are significant
  mutate(significant = !between(0, ci_lower, ci_upper),
         annotation = ifelse(significant, "*", ""))


# Add the correlations back to the data
trait_pop_str1 <- left_join(trait_pop_str, trait_pop_str_corr)

## Common plot modifier
g_mod <- list(
  geom_point(),
  geom_smooth(method = "lm", se = FALSE, color = "black"),
  geom_text(aes(x = Inf, y = -Inf, label = annotation), hjust = 1.5, vjust = 0.25, color = "black", size = 6),
  facet_grid(trait ~ PC, scales = "free"),
  scale_color_discrete(name = "Breeding\nProgram"),
  xlab("Eigenvalue"),
  theme_bw())


# Plot pop str vs mean
g_PC_v_mean <- trait_pop_str1 %>% 
  filter(measure == "g") %>%
  ggplot(aes(x = eigenvalue, y = value, col = program)) + 
  ylab("Genotypic Value") +
  g_mod

# Plot pop str vs linear stability
g_PC_v_b <- trait_pop_str1 %>% 
  filter(measure == "b") %>%
  ggplot(aes(x = eigenvalue, y = value, col = program)) + 
  ylab("Linear Stability") +
  g_mod

# Plot pop str vs nonlinear stability
g_PC_v_delta <- trait_pop_str1 %>%
  filter(measure == "log_delta") %>%
  ggplot(aes(x = eigenvalue, y = value, col = program)) + 
  ylab("Non-Linear Stability") +
  g_mod

# Save
save_file <- file.path(fig_dir, "population_structure_versus_mean.jpg")
ggsave(filename = save_file, plot = g_PC_v_mean, height = 4, width = 6)  

save_file <- file.path(fig_dir, "population_structure_versus_linear_stability.jpg")
ggsave(filename = save_file, plot = g_PC_v_b, height = 4, width = 6)  

save_file <- file.path(fig_dir, "population_structure_versus_nonlinear_stability.jpg")
ggsave(filename = save_file, plot = g_PC_v_delta, height = 4, width = 6)  




### Association Results - Genomewide Scan

## Load the results of the genomewide scan
load(file.path(result_dir, "S2MET_pheno_fw_mean_gwas_results.RData"))


## Common plot modifiers
# QQ plot
g_mod_qq <- list(
  geom_ribbon(aes(x = p_exp_neg_log10, ymin = ci_upper, ymax = ci_lower), fill = "grey75",
              inherit.aes = FALSE),
  geom_abline(slope = 1, intercept = 0),
  geom_point(),
  facet_grid(trait ~ .),
  ylab(expression(Observed~-log[10](italic(p)))),
  xlab(expression(Expected~-log[10](italic(p)))),
  theme_bw()
)

# Manhattan plot
g_mod_man <- list(
  geom_point(),
  geom_hline(yintercept = -log10(alpha), lty = 2),
  # geom_hline(aes(yintercept = neg_log10_fdr10, lty = "FDR 10%")),
  scale_color_manual(values = color, guide = FALSE),
  ylab(expression(-log[10](italic(q)))),
  xlab("Position (Mbp)"),
  theme_bw(),
  theme_manhattan(),
  theme(panel.border = element_blank())
)


#### Mean *per se* genomewide scan  ####

## QQ Plot
gwas_pheno_mean_fw_adj_qq <- gwas_pheno_mean_fw_adj %>% 
  arrange(trait, coef, pvalue) %>%
  group_by(trait, coef) %>%
  mutate(p_exp = ppoints(n = n()), # Generate p values under null of no associations
         p_exp_neg_log10 = -log10(p_exp), # Convert to -log10(p)
         # Add a confidence interval based on the beta distribution (assumes independence of tests)
         ci_lower = -log10(qbeta(p = (alpha / 2), shape1 = seq(n()), rev(seq(n())))),
         ci_upper = -log10(qbeta(p = 1 - (alpha / 2), shape1 = seq(n()), rev(seq(n()))))) %>%
  select(trait, coef, plot_coef, marker, p_exp_neg_log10, pvalue_neg_log10, ci_lower, ci_upper)

g_mean_qq <- gwas_pheno_mean_fw_adj_qq %>% 
  filter(coef == "g") %>%
  ggplot(aes(x = p_exp_neg_log10, y = pvalue_neg_log10)) + 
  labs(title = "Genotype Mean GWAS QQ Plot") +
  g_mod_qq

# Save
save_file <- file.path(fig_dir, "gwas_pheno_genotype_mean_qq.jpg")
ggsave(filename = save_file, plot = g_mean_qq, width = 4, height = 7, dpi = 1000)




## Manhattan plot
g_mean_manhattan <- gwas_pheno_mean_fw_adj %>%
  filter(coef == "g") %>%
  mutate(color = if_else(chrom %in% seq(1, 7, 2), "B", "G")) %>%
  ggplot(aes(x = pos / 1000000, y = qvalue_neg_log10, group = chrom, col = color)) + 
  facet_grid(trait ~ chrom, switch = "x", scales = "free", space = "free_x") +
  labs(title = expression('Genomewide Scan for Mean'~italic('Per Se'))) +
  g_mod_man

save_file <- file.path(fig_dir, "gwas_pheno_genotype_mean_manhattan.jpg")
ggsave(filename = save_file, plot = g_mean_manhattan, width = 9, height = 6, dpi = 1000)





#### Stability QTL
## QQ plot
g_stab_qq <- gwas_pheno_mean_fw_adj_qq %>% 
  ungroup() %>%
  filter(coef != "g") %>%
  ggplot(aes(x = p_exp_neg_log10, y = pvalue_neg_log10)) + 
  labs(title = "Phenotypic Stability GWAS QQ Plot") +
  g_mod_qq +
  facet_grid(trait ~ plot_coef)

save_file <- file.path(fig_dir, "gwas_pheno_stability_qq.jpg")
ggsave(filename = save_file, plot = g_stab_qq, width = 6, height = 7, dpi = 1000)


## Manhattan plot
g_stab_manhattan <- gwas_pheno_mean_fw_adj %>%
  # filter(model == "G", coef != "g") %>%
  filter(coef != "g") %>%
  mutate(color = if_else(chrom %in% seq(1, 7, 2), "B", "G")) %>%
  ggplot(aes(x = pos / 1000000, y = qvalue_neg_log10, group = chrom, col = color)) + 
  facet_grid(trait + plot_coef ~ chrom, switch = "x", scales = "free", space = "free_x") +
  labs(title = 'Genomewide Scan for Phenotypic Stability') +
  g_mod_man

save_file <- file.path(fig_dir, "gwas_pheno_stab_manhattan.jpg")
ggsave(filename = save_file, plot = g_stab_manhattan, width = 9, height = 9, dpi = 1000)





## Plot the results of the multi-locus mixed model

# Find the trait-markers that are not in the 'gwas_mlmm_Gmodel' df
gwas_pheno_mean_fw_adj_nonsig <- gwas_pheno_mean_fw_adj %>% 
  select(marker, trait, coef) %>% 
  setdiff(., select(ungroup(gwas_mlmm_Gmodel), marker, trait, coef)) %>%
  # left_join(., filter(gwas_pheno_mean_fw_adj, model == "G")) %>%
  left_join(., gwas_pheno_mean_fw_adj) %>%
  mutate(qvalue_neg_log10 = 0) %>% 
  select(trait, coef, marker, chrom, pos, beta, pvalue, qvalue, qvalue_neg_log10) %>% 
  rename(alpha = beta)

# Calculate the neg-log q values for the 'gwas_mlmm_Gmodel' df
gwas_mlmm_adj <- gwas_mlmm_Gmodel %>% 
  select(trait:pos, pvalue = pvalue, qvalue) %>%
  group_by(trait, coef) %>%
  mutate(qvalue_neg_log10 = -log10(qvalue)) %>%
  ungroup()

# Combine with the significant markers
gwas_pheno_mean_fw_adj_mlmm <- bind_rows(gwas_pheno_mean_fw_adj_nonsig, gwas_mlmm_adj) %>%
  mutate(color = if_else(chrom %in% seq(1, 7, 2), "B", "G"),
         plot_coef = str_replace_all(coef, coef_replace))


# Plot the results of trait mean per se
g_man_g <- gwas_pheno_mean_fw_adj_mlmm %>%
  filter(coef == "g") %>%
  mutate(chrom = as.factor(chrom)) %>% 
  ggplot(aes(x = pos / 1000000, y = qvalue_neg_log10, group = chrom, col = color)) + 
  facet_grid(trait ~ chrom, switch = "x", scales = "free", space = "free_x") +
  labs(title = expression('Multi-Locus Model for Mean'~italic('Per Se'))) +
  g_mod_man


# Save
save_file <- file.path(fig_dir, "gwas_pheno_genotype_mean_manhattan_mlmm.jpg")
ggsave(filename = save_file, plot = g_man_g, width = 9, height = 6, dpi = 1000)



# Plot the results of stability
g_man_stab <- gwas_pheno_mean_fw_adj_mlmm %>%
  filter(coef != "g") %>%
  ggplot(aes(x = pos / 1000000, y = qvalue_neg_log10, group = chrom, col = color)) + 
  facet_grid(trait + plot_coef ~ chrom, switch = "x", scales = "free", space = "free_x") +
  labs(title = 'Multi-Locus Model for Phenotypic Stability') +
  g_mod_man

# Save
save_file <- file.path(fig_dir, "gwas_pheno_fw_manhattan_mlmm.jpg")
ggsave(filename = save_file, plot = g_man_stab, width = 9, height = 9, dpi = 1000)







## Adjust the manhattan plots to show the markers that were declared significant 
## in the mlmm with different colors

# Get the SNPs not in the mlmm df
gwas_not_sig <- anti_join(gwas_pheno_mean_fw_adj, gwas_mlmm_adj, 
                                by = c("trait", "coef", "marker", "chrom", "pos")) %>%
  # Annotate the color for plotting
  mutate(color = if_else(chrom %in% seq(1, 7, by = 2), "B", "G"))

# Annotate the color from the mlmm
gwas_mlmm_sig <- gwas_mlmm_adj %>%
  mutate(color = "Bl")

# Combine
gwas_results_toplot <- bind_rows(gwas_mlmm_sig, gwas_pheno_not_sig) %>%
  select(trait:qvalue_neg_log10, color) %>%
  mutate(plot_coef = str_replace_all(coef, coef_replace), 
         plot_coef = parse_factor(plot_coef, levels = c("Genotype Mean", "Linear Stability", 
                                                        "Non-Linear Stability")))


## Manhattan plot for the scan and MLMM
## Iterate over each trait
for (tr in unique(gwas_results_toplot$trait)) {
  
  g_gwas_complete_plot <- gwas_results_toplot %>% 
    filter(trait == tr) %>%
    ggplot(aes(x = pos / 1000000, y = qvalue_neg_log10, group = chrom, color = color)) + 
    geom_point() + 
    g_mod_man +
    facet_grid(trait + plot_coef ~ chrom, switch = "x", scales = "free", space = "free_x")
  
  save_file <- file.path(fig_dir, str_c("gwas_manhattan_complete_", tr, ".jpg"))
  ggsave(filename = save_file, plot = g_gwas_complete_plot, width = 9, height = 6, dpi = 1000)
  
}



#### Summary of significant associations

## What are the numbers of significant marker-trait assocations for each trait and type?
gwas_sig_snp_summ <- gwas_pheno_mean_fw_adj %>% 
  group_by(trait, coef, chrom) %>% 
  summarize(n_sig_SNP = sum(qvalue <= alpha)) %>%
  mutate(type = "GWAS")

mlmm_sig_snp_summ <- gwas_mlmm_Gmodel %>% 
  select(trait:pos, qvalue) %>%
  ungroup() %>%
  complete(trait, coef, chrom) %>%
  group_by(trait, coef, chrom) %>% 
  summarize(n_sig_SNP = sum(qvalue <= alpha, na.rm = T)) %>%
  mutate(type = "MLMM")

# Combine with the multi-locus mixed model results
gwas_pheno_mean_fw_sig <- bind_rows(gwas_sig_snp_summ, mlmm_sig_snp_summ)

# Sumarize
gwas_pheno_mean_fw_sig %>% 
  group_by(trait, coef, type) %>% 
  summarize(n_sig_SNP = sum(n_sig_SNP)) %>% 
  spread(type, n_sig_SNP)



# Look at the significant markers and assess the effect size, minor allele frequency, etc.

# Calculate the allele frequency of the 1 allele
af1 <- {colMeans(M + 1) / 2} %>%
  data.frame(marker = names(.), af = ., row.names = NULL, stringsAsFactors = FALSE)

# Combine the minor allele frequency information
gwas_mlmm_marker_info <- gwas_mlmm_Gmodel %>% 
  ungroup() %>%
  left_join(., af1)



# Subset the marker matrix for these SNPs and count the number of lines
# from each program with each allele
gwas_mlmm_marker_prop <- gwas_mlmm_marker_info %>%
  group_by(marker) %>%
  do({
    df <- .
    
    # Extract the data on the marker
    marker_geno <- M[,df$marker, drop = FALSE] %>% 
      as.data.frame() %>% rownames_to_column("line_name") %>% 
      mutate(program = str_extract(line_name, "[A-Z][A-Z0-9]")) %>%
      select(line_name, program, marker = df$marker)
    
    # Round the marker genotype
    marker_geno_round <- marker_geno %>% 
      mutate(allele = case_when(marker <= -0.5 ~ -1,
                                marker >= 0.5 ~ 1,
                                TRUE ~ 0),
             allele = allele + 1)
    
    # Calculate the % of lines from each program that have the 1 allele
    marker_geno_prop <- marker_geno_round %>% 
      group_by(program) %>% 
      summarize(prop_1_allele = mean(allele) / 2) %>% 
      spread(program, prop_1_allele) 
    
    # Return the df
    bind_cols(df, marker_geno_prop)
    
  })

gwas_mlmm_marker <- gwas_mlmm_marker_prop %>%
  ungroup() %>%
  left_join(., snp_info) %>%
  select(trait, coef, marker, chrom, pos, cM = cM_pos, alpha, qvalue, R_sqr_snp, af, AB:WA) %>%
  mutate(coef = str_replace_all(coef, coef_replace),
         MAF = pmin(af, 1 - af)) %>%
  select(trait:cM, MAF, alpha:R_sqr_snp, AB:WA) %>%
  arrange(trait, coef, chrom, pos)





### Test for pleiotropy

# Load the results
load(file.path(result_dir, "/S2MET_pheno_fw_gwas_pleiotropy_results.RData"))

## Tidy up the data
gwas_pheno_mean_fw_plei <- gwas_pheno_mean_fw_plei_out %>%
  list(., names(.)) %>%
  pmap_df(~mutate(.x, trait = .y)) %>%
  separate(trait, c("trait", "stab_coef"), sep = "\\.") %>%
  mutate(term = if_else(term == "g", term, stab_coef)) %>%
  select(marker, trait, term, stab_coef, beta, se)

## The test for pleiotropy is a union-intersection test
## For each marker, there are two betas: beta_g and beta_stability,
## where beta_g is the marker effect for the genotype mean and beta_stability
## is the marker effect for phenotypic stability.
## 
## H0: beta_g == 0 or beta_stability == 0
## HA: beta_g != 0 and beta_stability != 0
## 
## We can use the maximum p-value as the test - if this p-value is less
## than a threshold alpha, the null hypothesis of no pleiotropy can be rejected.
## 

# Calculate test statistics and p-values
gwas_pheno_mean_fw_plei_test <- gwas_pheno_mean_fw_plei %>% 
  mutate(statistic = (beta^2) / (se^2), 
         pvalue = pchisq(q = statistic, df = 1, lower.tail = FALSE))

# Correct the p-values for association test (i.e. mean per se and stability)
gwas_pheno_mean_fw_plei_adj <- gwas_pheno_mean_fw_plei_test %>% 
  group_by(trait, term, stab_coef) %>% 
  mutate(pvalue_adj = p.adjust(pvalue, "fdr"),
         qvalue = qvalue(pvalue)$qvalue) %>%
  group_by(marker, trait, stab_coef) %>%
  summarize_at(vars(pvalue:qvalue), max) %>%
  ungroup()

## None of the markers were above the alpha threshold. Instead, use an empirical 
## threshold of the top 0.1% of markers for each test
gwas_pheno_mean_fw_plei_adj_sig <- gwas_pheno_mean_fw_plei_adj %>% 
  group_by(trait, stab_coef) %>% 
  # Use less-than for instances like plant height non-linear stability, where
  # no q-value is below 1.
  mutate(significant = qvalue < quantile(qvalue, 0.001)) 
  

# Get the allele effects
gwas_pheno_mean_fw_plei_adj_sig_effects <- gwas_pheno_mean_fw_plei_adj_sig %>%
  ungroup() %>%
  filter(significant) %>%
  left_join(., gwas_pheno_mean_fw_plei_test, by = c("marker", "trait", "stab_coef")) %>%
  left_join(., snp_info) %>%
  select(marker, chrom, pos, trait, stability = stab_coef, term, effect = beta) %>%
  mutate(stability = str_replace_all(stability, coef_replace),
         term = if_else(term == "g", "genotype_mean_marker_effect", "stability_marker_effect")) %>% 
  spread(term, effect) %>%
  arrange(trait, stability, chrom, pos)




# Combine the data add SNP information
gwas_pheno_mean_fw_plei_toplot <- gwas_pheno_mean_fw_plei_adj_sig %>%
  ungroup() %>%
  full_join(., snp_info, by = "marker") %>%
  select(marker, chrom:cM_pos, names(.)) %>%
  arrange(trait, stab_coef, chrom, pos) %>%
  mutate(
    stab_coef = str_replace_all(stab_coef, coef_replace),
    stab_coef = parse_factor(stab_coef, levels = c("Linear Stability", "Non-Linear Stability")),
    color = if_else(chrom %in% seq(1, 7, 2), "B", "G"),
    color = if_else(significant, "Or", color))

# Labeller function
label_coef <- function(x) str_c("Genotype Mean\n", x)


## Manhattan plot for the pleitropic scan
## Iterate over each trait
for (tr in unique(gwas_pheno_mean_fw_plei_toplot$trait)) {
  
  g_gwas_plei <- gwas_pheno_mean_fw_plei_toplot %>%
    filter(trait == tr) %>%
    ggplot(aes(x = pos / 1000000, y = -log10(qvalue), color = color)) + 
    geom_point() + 
    g_mod_man +
    facet_grid(trait + stab_coef ~ chrom, scales = "free", space = "free_x", switch = "x",
               labeller = labeller(stab_coef = label_coef))
  
  save_file <- file.path(fig_dir, str_c("gwas_plei_manhattan_", tr, ".jpg"))
  ggsave(filename = save_file, plot = g_gwas_plei, width = 9, height = 6)
  
}







#### Resampling the Number of Environments

## Resample 20, 40, 60, or 80% of environments 250 times each, calculate the stability 
## estimates, then performed the mapping.

# Read in the results
load(file.path(result_dir, "S2MET_pheno_fw_gwas_resample_results.RData"))

# Bind rows
resample_gwas_sig_out1 <- resample_gwas_sig_out %>%
  bind_rows()

## Filter the original GWAS results for the stability QTL and the G model
gwas_sig_fw <- gwas_pheno_mean_fw_adj %>% 
  filter(coef != "g", 
         q_value <= alpha)

# Iterate over the resampling results and find the signficant loci at alpha
resample_gwas_sig_fw <- resample_gwas_sig_out1 %>%
  select(p, iter, sig_out) %>% 
  unnest() %>%
  # SPlit the trait into trait and coefficient
  separate(col = trait, into = c("trait", "coef"), sep = "_", extra = "merge") %>% 
  filter(q_value <= alpha)

## First, for each proportion of environments and for each level of stability,
## how many significant SNPs were detected?
resample_gwas_sig_fw_count <- resample_gwas_sig_fw %>% 
  group_by(p, iter, trait, coef) %>% 
  summarize(n_SNP = n()) %>%
  ungroup() %>%
  complete(p, iter, trait, coef, fill = list(n_SNP = 0)) %>%
  group_by(p, trait, coef) %>%
  mutate_at(vars(n_SNP), funs(mean, sd)) %>%
  ungroup()

## Plot this
g_resample_n_SNPs <- resample_gwas_sig_fw_count %>% 
  ggplot(aes(x = p, y = mean, col = trait, group = trait)) + 
  geom_point() + 
  geom_line() + 
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.25) + 
  ylab("Number of Significant SNPs") +
  xlab("Proportion of Environments Sampled") +
  facet_grid(~ coef) + 
  theme_bw()




# Subset the original gwas results
gwas_sig_fw_intersect <- gwas_sig_fw %>%
  select(marker, trait, coef)

# For each p, iter, trait, and coef, how often do we detect new SNPs or the same
# SNPs as in the full analysis?
resample_gwas_sig_fw_intersect <- resample_gwas_sig_fw %>% 
  select(p, iter, marker, trait, coef) %>%
  group_by(p, iter, trait, coef) %>%
  do({
    df <- .
    
    # Trait and coef
    tr <- unique(df$trait)
    co <- unique(df$coef)
    
    # Filter the original data
    gwas_sig_fw_intersect1 <- gwas_sig_fw_intersect %>%
      filter(trait == tr, coef == co)
    
    # Intersect
    marker_intersect <- intersect(df$marker, gwas_sig_fw_intersect1$marker)
    marker_new <- setdiff(df$marker, gwas_sig_fw_intersect1$marker)
    
    # Return a data.frame of the intersected markers, the unique markers, and the counts of each
    data_frame(hit_SNPs = list(marker_intersect), 
               new_SNPs = list(marker_new))  })

resample_gwas_sig_fw_intersect1 <- resample_gwas_sig_fw_intersect %>% 
  ungroup() %>% 
  gather(SNP_type, SNPs, -p:-coef) %>% 
  mutate(n_SNPs = map_dbl(SNPs, length)) %>% 
  group_by(p, trait, coef, SNP_type) %>% 
  summarize_at(vars(n_SNPs), funs(mean, sd)) %>%
  ungroup()

## Remember, the above measures the number of intersected or new SNPs conditional
## on actually detected significant associations

# Number of iterations
n_iter <- n_distinct(resample_gwas_sig_fw$iter)

# Create a data.frame to plot
resample_gwas_sig_fw_intersect_toplot <- resample_gwas_sig_fw_intersect1 %>% 
  mutate(coef = str_replace_all(coef, coef_replace), 
         SNP_type = str_replace_all(SNP_type, c("hit_SNPs" = "Recovered SNPs", "new_SNPs" = "New SNPs")))


# Plot
g_resample_intersect <- resample_gwas_sig_fw_intersect_toplot %>% 
  ggplot(aes(x = p, y = mean, fill = SNP_type, group = SNP_type)) + 
  geom_col(position = "dodge") + 
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.25, 
                position = position_dodge(0.9)) + 
  ylab("Number of SNPs") +
  xlab("Proportion of Environments Sampled") +
  scale_fill_discrete(name = "Association\nType") + 
  facet_grid(trait ~ coef) + 
  theme_bw() +
  theme(legend.position = "bottom")

## Count up the intersected and new SNPs and rank by number of times
## the SNP was detected

resample_snp_detect_count <- resample_gwas_sig_fw_intersect %>% 
  ungroup() %>% 
  gather(SNP_type, marker, -p:-coef) %>%
  unnest() %>% 
  group_by(p, trait, coef, SNP_type, marker) %>% 
  summarize(times_detected = n(), prop_detected = times_detected / n_iter) %>%
  ungroup() %>%
  left_join(., select(S2TP_imputed_multi_genos_hmp, marker = rs, chrom, pos), by = "marker") %>%
  # Assign proportion bins 
  mutate(prop_group = cut(x = prop_detected, breaks = seq(0, 1, 0.2), include.lowest = TRUE, right = TRUE))

# Create a data.frame for plotting
resample_snp_detect_count_toplot <- resample_snp_detect_count %>%
  mutate_at(vars(p, trait, coef, SNP_type), funs(as.factor)) %>% 
  complete(p, trait, coef, SNP_type, fill = list(prop_detected = 0)) %>%
  mutate(coef = str_replace_all(coef, coef_replace), 
         SNP_type = str_replace_all(SNP_type, c("hit_SNPs" = "Recovered SNPs", "new_SNPs" = "New SNPs")))

# Are we more likely to discover new SNPs or those we discovered before?
g_resample_snp_detect_count <- resample_snp_detect_count_toplot %>%
  ggplot(aes(x = p, y = prop_detected, color = SNP_type)) +
  geom_boxplot(position = "dodge", width = 0.5) +
  scale_color_discrete(drop = FALSE) +
  facet_grid(trait ~ coef) + 
  ylab("Probability of Detecting the SNP") + 
  xlab("Proportion of Environments Sampled") +
  scale_color_discrete(name = "Association\nType") + 
  theme_bw() +
  theme(legend.position = "bottom")


## Combine the plots
# g_resample_combined <- g_resample_intersect + g_resample_snp_detect_count
g_resample_combined <- plot_grid(g_resample_intersect, g_resample_snp_detect_count, 
                                 labels = c("A", "B"))


save_file <- file.path(fig_dir, "resample_gwas_detection_results.jpg")
ggsave(filename = save_file, plot = g_resample_combined, width = 8, height = 5)





# Sort and display
resample_snp_detect_count_max <- resample_snp_detect_count %>%
  filter(times_detected == max(times_detected))


# Common plot modifier
g_mod <- list(geom_point(),
              geom_segment(data = barley_lengths, aes(x = 0, y = 0, xend = 0, yend = length / 1000000), 
                           inherit.aes = FALSE, size = 2, lineend = "round", col = "grey70") ,
              ylab("Position (Mbp)") ,
              xlab("Chromosome"),
              xlim(c(-2, 2)) , 
              facet_grid(~ chrom, switch = "x") , 
              scale_shape_manual(values = c("Recovered SNPs" = 1, "New SNPs" = 0), name = "Association\nType"), 
              scale_size_discrete(name = "Probability of\nDetecting SNP"),
              scale_color_discrete(name = "Stability\nType"),
              scale_y_reverse() , 
              theme_bw() , 
              theme(panel.grid = element_blank(), 
                    panel.border = element_blank(), 
                    # panel.spacing.x = unit(1, units = "cm"),
                    axis.text.x = element_blank(), 
                    axis.title.x = element_blank(), 
                    axis.ticks.x = element_blank(),
                    strip.background = element_blank()) )


# Create a data.frame to plot
resample_snp_detect_count_toplot <- resample_snp_detect_count %>%
  mutate(p = parse_number(p),
         x = if_else(SNP_type == "hit_SNPs", -p * 2, p * 2),
         coef = str_replace_all(coef, coef_replace),
         SNP_type = str_replace_all(SNP_type, c("hit_SNPs" = "Recovered SNPs", "new_SNPs" = "New SNPs")))



# Create a 'heatmap' for each trait
g_resample_hd <- resample_snp_detect_count_toplot %>% 
  filter(trait == "HeadingDate") %>% 
  ggplot(aes(y = pos / 1000000, x = x, size = prop_group, col = coef, shape = SNP_type)) + 
  g_mod +
  labs(title = "Heading Date")

g_resample_ph <- resample_snp_detect_count_toplot %>% 
  filter(trait == "PlantHeight") %>% 
  ggplot(aes(y = pos / 1000000, x = x, size = prop_group, col = coef, shape = SNP_type)) + 
  g_mod +
  labs(title = "Plant Height")

g_resample_gy <- resample_snp_detect_count_toplot %>% 
  filter(trait == "GrainYield") %>% 
  ggplot(aes(y = pos / 1000000, x = x, size = prop_group, col = coef, shape = SNP_type)) + 
  g_mod +
  labs(title = "Grain Yield")

## Plot together
save_file <- file.path(fig_dir, "resample_count_hd.jpg")
ggsave(filename = save_file, plot = g_resample_hd, width = 10, height = 7, dpi = 1000)

save_file <- file.path(fig_dir, "resample_count_ph.jpg")
ggsave(filename = save_file, plot = g_resample_ph, width = 10, height = 7, dpi = 1000)

save_file <- file.path(fig_dir, "resample_count_gy.jpg")
ggsave(filename = save_file, plot = g_resample_gy, width = 10, height = 7, dpi = 1000)







## Save data for further analysis
save_file <- file.path(result_dir, "gwas_adjusted_significant_results.RData")
save("gwas_pheno_mean_fw_adj", "gwas_mlmm_marker_info", "gwas_pheno_mean_fw_plei_toplot", 
     file = save_file)





