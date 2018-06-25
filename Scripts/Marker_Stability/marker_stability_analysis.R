## S2MET Mapping
## Calculate marker effect reaction norms
## 
## 
## Author: Jeff Neyhart
## Last updated: June 14, 2018
## 

# Extra libraries
library(cowplot)
library(lme4qtl)

# Load the source
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))
 
# Load the stability results
load(file.path(result_dir, "pheno_mean_fw_results.RData"))
# Load the marker effect by environment results
load(file.path(result_dir, "marker_by_env_effects.RData"))
# Load the GWAS results
load(file.path(result_dir, "pheno_fw_mean_gwas_results.RData"))


# Get the marker information
snp_info <- S2TP_imputed_multi_genos_hmp %>%
  select(marker = rs, chrom, pos, cM_pos) %>%
  # Correct the cM position for BOPA snps
  mutate(cM_pos = if_else(str_detect(marker, "^S"), cM_pos, cM_pos / 1000))


M <- S2TP_imputed_multi_genos_mat


# GWAS model
model_use <- "QG"
# Significance
alpha <- 0.05



# ## Manage the mxe data and convert to a tidy data.frame
# marker_by_env_effects_df <- marker_by_env_effects %>% 
#   map(~do.call("cbind", .)) %>%
#   map(., ~t(.) %>% as.data.frame() %>% rownames_to_column("marker") %>%
#         rename_at(vars(-marker), ~str_replace(string = ., "environment", "")) %>%
#         gather(environment, effect, -marker) %>% as_data_frame() )
#   
# mxe_df <- marker_by_env_effects_df %>%
#   list(., names(.)) %>% 
#   pmap_df(~mutate(.x, trait = .y)) %>%
#   select(trait, environment, names(.))

mxe_df <- marker_by_env_effects %>%
  unnest(marker_effect) %>%
  ungroup()

# Combine with the environmental mean
mxe_df1 <- mxe_df %>%
  left_join(., distinct(S2_MET_pheno_mean_fw, trait, environment, h),
            by = c("environment", "trait"))

# Calculate the marker effect stability
# Regress the marker effects in each environment on the mean in that environment
marker_effect_stab <- mxe_df1 %>%
  group_by(trait, marker) %>%
  do({
    df <- .
    # Fit the model
    fit <- lm(effect ~ h, df)
    
    # Output a data.frame with the intercept, slope, MSE, and fitted marker effects
    mutate(df, a = coef(fit)[1], c = coef(fit)[2], delta = mean(resid(fit)^2), 
           fitted_effects = fitted(fit)) 
    
  })


# Add the snp information
marker_mean_fw <- marker_effect_stab %>%
  left_join(., snp_info, by = "marker") %>%
  select(trait, marker, chrom:cM_pos, trait, names(.)) %>%
  ungroup()

# Save the data
save("marker_mean_fw", file = file.path(result_dir, "marker_mean_fw_results.RData"))



# ## Are polynomial regressions of marker effects better fits?
# marker_effect_stab_fits <- mxe_df1 %>%
#   group_by(trait, marker) %>%
#   do({
#     df <- .
#     
#     # Fit the full model with polynomial regression
#     fit_full <- lm(effect ~ 1 + poly(h, 2), df)
#     # Fit the linear regression
#     fit_lin <- lm(effect ~ 1 + h, df)
#     # Fit the base
#     fit_base <- lm(effect ~ 1, df)
#     
#     aic_compare <- AIC(fit_base, fit_lin, fit_full)
#     best <- subset(aic_compare, AIC == min(AIC))
#     
#     # Return the best model
#     data_frame(fit_type = row.names(best))
#   })
# 
# # Summarize by trait
# marker_effect_stab_fits %>% 
#   group_by(trait, fit_type) %>% 
#   summarize(n = n()) %>% 
#   spread(fit_type, n) %>%
#   mutate_at(vars(fit_base:fit_lin), ~ . / ncol(M))
# 





###### Below here, data from above can be loaded

load(file.path(result_dir, "marker_mean_fw_results.RData"))

## Calculate the mean marker effect of each marker
mean_marker_effect <- marker_mean_fw %>% 
  group_by(marker, trait) %>% 
  mutate(mean_marker_effect = mean(effect)) %>%
  ungroup()

# Plot
mean_marker_effect %>% 
  distinct(marker, trait, mean_marker_effect) %>% 
  qplot(x = mean_marker_effect, data = .) + 
  facet_wrap(~ trait, ncol = 2, scales = "free_x") + 
  theme_bw()

# Scan
mean_marker_effect %>% 
  distinct(marker, chrom, pos, trait, mean_marker_effect) %>%
  ggplot(aes(x = pos / 1000000, y = mean_marker_effect)) + 
  geom_point() + 
  ylab("Marker Effect Plasticity Estimate") +
  xlab("Position (Mbp)") +
  facet_grid(trait ~ chrom, scales = "free", space = "free_x", switch = "both",
             labeller = labeller(chrom = function(x) str_c("Chr. ", x)))   +
  theme_manhattan()




## Plot the distribution of stability estimates

## Transform
marker_mean_fw_trans <- marker_mean_fw %>% 
  distinct(marker, chrom, pos, trait, c, delta) %>%
  mutate(log_delta = log(delta), c = c + 1)


## Plot the stability terms individually, then combine
# First define a common list of ggplot modifiers
g_mod <- list(
  geom_density(aes(fill = "blue")),
  xlab("Estimate"),
  theme_bw() +
  theme(axis.title.y = element_blank()) )

# Just plot linear stability
g_marker_fw_dens_c <- marker_mean_fw_trans %>%
  ggplot(aes(x = c)) + 
  labs(title = "Linear Stability") +
  facet_wrap( ~ trait, ncol = 1) +
  scale_fill_brewer(palette = "Set2", guide = FALSE) +
  g_mod


# Just plot non-linear stability
g_marker_fw_dens_delta <- marker_mean_fw_trans %>%
  ggplot(aes(x =  log_delta)) + 
  labs(title = "Non-Linear Stability") +
  scale_fill_brewer(palette = "Set2", guide = FALSE) +
  facet_wrap( ~ trait, ncol = 1, scales = "free_x") +
  g_mod

# Add the plots together
g_marker_fw_dist <- plot_grid(g_marker_fw_dens_c, g_marker_fw_dens_delta, ncol = 2)

ggsave(filename = "marker_stability_estimate_distriubtions.jpg", plot = g_marker_fw_dist, 
       height = 6, width = 5, dpi = 1000, path = fig_dir)


# Plot for all markers the relationship between linear and nonlinear stability
g_linear_vs_nonlinear <- marker_mean_fw_trans %>% 
  qplot(x = c, y = log_delta, data = .) + 
  facet_wrap(~trait, ncol = 2, scales = "free_y") +
  ylab("Non-Linear Stability") +
  xlab("Linear Stability") +
  theme_bw()

ggsave(filename = "marker_linear_vs_nonlinear_stability.jpg", plot = g_linear_vs_nonlinear, 
       height = 6, width = 5, dpi = 1000, path = fig_dir)

# Plot for all markers the relationship between the slope and the intercept
marker_mean_fw %>% 
  distinct(trait, marker, a, c) %>%
  qplot(x = a, y = c, data = .) +
  facet_wrap(~trait, ncol = 2, scales = "free_x") +
  theme_bw()

## Correlate
marker_mean_fw %>% 
  distinct(trait, marker, a, c)  %>% 
  group_by(trait) %>% 
  summarize(cor = cor(a, c))



## Assess what markers change sign across environments`
marker_effect_sign <- marker_mean_fw %>% 
  mutate(effect_sign = sign(effect)) %>% 
  group_by(trait, marker, effect_sign) %>% 
  summarize(n = n()) %>% 
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  complete(trait, marker, effect_sign, fill = list(n = 0, prop = 0))

## Proportion of markers that change sign
marker_effect_sign %>% 
  filter(effect_sign != -1) %>% 
  group_by(trait) %>% 
  summarize(n_change_sign = sum(prop != 1 & prop != 0), 
            prop_change_sign = n_change_sign / n())

## Plot a histogram of the proportion of positive effects
marker_effect_sign %>% 
  filter(effect_sign != -1) %>% 
  qplot(x = prop, data = ., geom = "density") + 
  facet_wrap(~trait, ncol = 2) +
  theme_bw()



## Take the significant markers from GWAS and look at their effect reaction norms
gwas_mlmm_model_adj <- subset(gwas_mlmm_model, model == model_use) %>%
  mutate(mlmm_out = map(mlmm_out, "fit_out_reduced") %>% 
           map(~mutate(subset(.$summary, term != "Q"), snp_r_squared = .$r_squared$snp_r_squared, 
                       all_snp_r_squared = .$r_squared$fixed_r_squared["all_snps"]))) %>%
  unnest(mlmm_out)


# Plot for all markers the relationship between linear and nonlinear stability
g_linear_vs_nonlinear <- marker_mean_fw_trans %>% 
  qplot(x = c, y = log_delta, data = .) + 
  facet_wrap(~trait, ncol = 2, scales = "free_y") +
  ylab("Non-Linear Stability") +
  xlab("Linear Stability") +
  theme_bw()

ggsave(filename = "marker_linear_vs_nonlinear_stability.jpg", plot = g_linear_vs_nonlinear, 
       height = 6, width = 5, dpi = 1000, path = fig_dir)

# Plot for all markers the relationship between the slope and the intercept
marker_mean_fw %>% 
  distinct(trait, marker, a, c) %>%
  qplot(x = a, y = c, data = .) +
  facet_wrap(~trait, ncol = 2, scales = "free_x") +
  theme_bw()

## Correlate
marker_mean_fw %>% 
  distinct(trait, marker, a, c)  %>% 
  group_by(trait) %>% 
  summarize(cor = cor(a, c))



## Assess what markers change sign across environments`
marker_effect_sign <- marker_mean_fw %>% 
  mutate(effect_sign = sign(effect)) %>% 
  group_by(trait, marker, effect_sign) %>% 
  summarize(n = n()) %>% 
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  complete(trait, marker, effect_sign, fill = list(n = 0, prop = 0))

## Proportion of markers that change sign
marker_effect_sign %>% 
  filter(effect_sign != -1) %>% 
  group_by(trait) %>% 
  summarize(n_change_sign = sum(prop != 1 & prop != 0), 
            prop_change_sign = n_change_sign / n())

## Plot a histogram of the proportion of positive effects
marker_effect_sign %>% 
  filter(effect_sign != -1) %>% 
  qplot(x = prop, data = ., geom = "density") + 
  facet_wrap(~trait, ncol = 2) +
  theme_bw()

marker_mean_fw_gwas_mlmm <- select(marker_mean_fw, trait, marker, environment, effect, h) %>%
  left_join(select(gwas_mlmm_model_adj, trait, coef, marker = term), .)

## Take the significant markers from GWAS and look at their effect reaction norms
gwas_mlmm_model_adj <- subset(gwas_mlmm_model, model == model_use) %>%
  mutate(mlmm_out = map(mlmm_out, "fit_out_reduced") %>% 
           map(~mutate(subset(.$summary, term != "Q"), snp_r_squared = .$r_squared$snp_r_squared, 
                       all_snp_r_squared = .$r_squared$fixed_r_squared["all_snps"]))) %>%
  unnest(mlmm_out)

marker_mean_fw_gwas_mlmm <- select(marker_mean_fw, trait, marker, environment, effect, h) %>%
  left_join(select(gwas_mlmm_model_adj, trait, coef, marker = term), .)


# Plot the reaction norms
g_marker_fw_gwas <- marker_mean_fw_gwas_mlmm %>% 
  ggplot(aes(x = h, y = effect, color = environment, group = 1)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(~ trait + coef + marker, ncol = 3, scales = "free") +
  theme_bw()

ggsave(filename = "marker_fw_gwas_rn.jpg", plot = g_marker_fw_gwas, path = fig_dir,
       height = 10, width = 6, dpi = 1000)


# First group markers into those that are highly stable, highly sensitive, or neither using empirical thresholds

# Tidy up
marker_mean_fw_tidy <- marker_mean_fw_trans %>% 
  select(-delta) %>% 
  gather(coef, estimate, -trait:-pos)
  

# What should be the cutoff level?
alpha <- 0.05

# For each trait, calculate empirical thresholds for significance
marker_mean_fw_sig <- marker_mean_fw_tidy %>%
  filter(coef == "c") %>% 
  group_by(trait) %>% 
  mutate(lower_perc = quantile(estimate, alpha / 2), 
         upper_perc = quantile(estimate, 1 - (alpha / 2))) %>%
  ungroup() %>%
  mutate(significance = case_when(estimate >= upper_perc ~ "Plastic",
                                  estimate <= lower_perc ~ "Plastic",
                                  TRUE ~ "Stable"),
         marker_type = if_else(str_detect(marker, "^S"), "GBS", "BOPA"))

# # Calculate a table of markers and sensitive/stable/notsignificant
# mar_stab_table <- S2_MET_marker_eff_pheno_fw_sig %>%
#   group_by(trait, marker_type, significance) %>%
#   summarize(n = n()) %>%
#   mutate(prop = n / sum(n))
# 
# mar_stab_table %>% 
#   select(-n) %>% 
#   spread(significance, prop)


## Table of traits and marker significance groups
xtabs(~ trait + significance, marker_mean_fw_sig)




## Plot the stability estimates across the genome with the empirical threshold
# Color theme
colors <- c("black", umn_palette(n = 4)[-c(1:3)]) %>%
  set_names(c("Stable", "Plastic"))


g_mar_stab <- marker_mean_fw_sig %>%   
  ggplot(aes(x = pos / 1000000, y = estimate, col = significance)) + 
  geom_point() + 
  geom_hline(aes(yintercept = lower_perc), lty = 2) +
  geom_hline(aes(yintercept = upper_perc), lty = 2) +
  scale_color_manual(values = colors, name = "Marker Type") +
  ylab("Marker Effect Plasticity Estimate") +
  xlab("Position (Mbp)") +
  facet_grid(trait ~ chrom, scales = "free", space = "free_x", switch = "both",
             labeller = labeller(chrom = function(x) str_c("Chr. ", x)))   +
  theme_manhattan()

# Save
ggsave(filename = "marker_stability_genome.jpg", plot = g_mar_stab, path = fig_dir,
       height = 6.5, width = 9, dpi = 1000)


# ## Poster version
# g_mar_stab_poster <- g_mar_stab +
#   geom_point(size = 2) + 
#   labs(caption = "95% empirical threshold given by dotted line.") +
#   theme_poster() + 
#   theme(legend.position = "bottom",
#         panel.grid = element_blank(),
#         axis.text.x = element_text(angle = 45, hjust = 1),
#         legend.text = element_text(size = 14),
#         panel.spacing.x = unit(x = 0, units = "cm"),
#         panel.spacing.y = unit(x = 1, units = "lines"),
#         panel.border = element_rect(color = "grey75"),
#         plot.caption = element_text(size = 10))
# 
# ggsave(filename = "marker_stability_genome_poster.jpg", plot = g_mar_stab_poster, path = fig_dir,
#        height = 10, width = 16, dpi = 1000)





#### Stable versus plastic markers

## Gene Annotation

# Look at gene annotation with relation to the different groups of SNPs. 
# Determine the proportion of SNPs in genes, gene-proximal, or intergenic


## Load the barley annotation
ann_dir <- "C:/Users/Jeff/GoogleDrive/BarleyLab/Projects/Genomics/Annotation"
load(file.path(ann_dir, "barley_genomic_ranges.RData"))


### First plot gene annotation density across the genome
barley_genes <- barley_grange_list$genes %>%
  subset(seqnames != "chrUn")
# Set the lengths of the seqnames to the lengths of the chromosomes
barley_seqlengths <- setNames(object = barley_lengths$length, 
                              nm = str_c("chr", barley_lengths$chrom, "H"))


seqlevels(barley_genes) <- head(seqlevels(barley_genes), -1)
seqlengths(barley_genes) <- barley_seqlengths

# Set a window size and step size
window_size <- step_size <- 10000000

## Split the annotation by chromosome
barley_genes_split <- split(x = barley_genes, f = seqnames(barley_genes)) %>% 
  as.list()

## Split the seqlengths by chromosome
barley_seqlengths_split <- list(names(barley_genes_split), seqlengths(barley_genes)) %>% 
  pmap(~tileGenome(seqlengths = setNames(.y, .x), tilewidth = window_size))

## For each chromosome, subset by overlap with the annotation list
barley_genes_perchrom <- list(barley_seqlengths_split, barley_genes_split) %>%
  pmap_df(~{
    
    # Find the overlaps between the genes and the sequences
    overlaps <- findOverlaps(subject = .y, query = .x) %>%
      as.data.frame()
    
    barley_seqlengths_df <- .x %>% 
      as.data.frame() %>% 
      select(group, seqnames:width)
    
    ## Merge and return
    left_join(overlaps, barley_seqlengths_df, by = c("queryHits" = "group")) %>%
      group_by(queryHits) %>% 
      do({ 
        df <- .
        df %>% 
          distinct(queryHits, seqnames, start, end, width) %>% 
          mutate(n_genes = n_distinct(.y[df$subjectHits, ]$gene_id))
      }) %>%
      ungroup()
  })
  
## Plot
g_gene_dens <- barley_genes_perchrom %>% 
  mutate(pos = (start + end) / 2,
         chrom = parse_number(seqnames)) %>%
  ggplot(aes(x = pos / 1e6, y = n_genes)) + 
  geom_smooth(method = "loess", span = 0.25) + 
  ylab("Number\nof Genes") +
  xlab("Position (Mbp)") +
  facet_grid(~ chrom, scales = "free_x", switch = "x", space = "free_x") +
  theme_manhattan()

## Combine the stability plot above with the gene density plot
g_plast_gene_dens <- plot_grid(g_mar_stab + theme(axis.title.x = element_blank()), 
                               g_gene_dens, ncol = 1, 
                               rel_heights = c(0.80, 0.20), align = "hv", axis = "lr")

## Save
ggsave(filename = "marker_stability_and_genes.jpg", plot = g_plast_gene_dens, path = fig_dir,
       height = 8, width = 10, dpi = 1000)









# Convert the SNPs in this experiment to a GRanges object
snp_info_adj <- S2TP_imputed_multi_genos_hmp %>%
  select(rs:cM_pos) %>%
  separate(alleles, c("ref", "alt"), "/") %>%
  mutate(chrom = str_c("chr", chrom, "H"))

snp_info_grange <- GRanges(seqnames = snp_info_adj$chrom,
                           ranges = IRanges(snp_info_adj$pos, snp_info_adj$pos),
                           ref = snp_info_adj$ref, alt = snp_info_adj$alt,
                           cM_pos = snp_info_adj$cM_pos, rs = snp_info_adj$rs)


# Find the distance to the nearest gene for each SNP
# Convert to a data.frame and add the SNP info
snp_info_distance_to_gene <- distanceToNearest(x = snp_info_grange, subject = barley_grange_list$genes) %>%
  as.data.frame() %>%
  bind_cols(snp_info_adj, .) %>%
  mutate(gene_id = barley_grange_list$genes[subjectHits]$gene_id) %>% # Add the gene id
  select(-queryHits, -subjectHits)

# Assign a cutoff for proximal
# The value below is in base-pairs
proximity_cutoff <- 10000

## Assign genic, proximal, or non-genic 
snp_info_distance_to_gene_ann <- snp_info_distance_to_gene %>%
  mutate(class = case_when(distance == 0 ~ "Genic",
                           between(distance, 0, proximity_cutoff) ~ "Proximal",
                           distance > proximity_cutoff ~ "Non-genic",
                           TRUE ~ as.character(NA)),
         class = parse_factor(class, levels = c("Genic", "Proximal", "Non-genic")),
         marker_type = if_else(str_detect(rs, "^[1-2]{2}"), "BOPA", "GBS")) %>%
  dplyr::rename(marker = rs)

# Find the proportion of SNPs that are in each group. This forms the null proportion
null_snp_prop <- snp_info_distance_to_gene_ann %>% 
  group_by(class) %>% 
  summarize(nSNP = n()) %>%
  mutate(prop = nSNP / sum(nSNP))

# Now use the marker stability estimates to group SNPs
snp_info_distance_to_gene_ann1 <- snp_info_distance_to_gene_ann %>%
  full_join(., select(marker_mean_fw_sig, marker, trait, significance),
            by = "marker")

# Group by the marker stability classes and calculate proportions
# Then add the null proportion
marker_stability_snp_prop <- snp_info_distance_to_gene_ann1 %>% 
  group_by(trait, significance, class) %>% 
  summarize(nSNP = n()) %>% 
  mutate(prop = nSNP / sum(nSNP),
         significance_SNP = sum(nSNP)) %>%
  left_join(., select(null_snp_prop, class, null_prop = prop), by = "class")
  

# Table
marker_stability_snp_prop %>%
  select(-nSNP, -null_prop) %>% 
  spread(class, prop)


## Set the seed for the binomial test
set.seed(512)

## Perform an binomial exact test
## H0: the proportion of SNPs in each significance group in each class is
## not different than the proportion across all SNPs
marker_stability_snp_prop_test <- marker_stability_snp_prop %>% 
  group_by(trait, significance, class) %>%
  mutate(binom_test = list({binom.test(x = nSNP, n = significance_SNP, p = null_prop)})) %>%
  ungroup() %>% 
  ## Extract p-value and confidence intervals
  mutate(out = map(binom_test, ~{ data.frame(p_value = .$p.value, estimate = .$estimate, 
                                             lower = .$conf.int[1], upper = .$conf.int[2]) }))
  
## Adjust the pvalues
marker_stability_snp_prop_adj <- marker_stability_snp_prop_test %>%
  unnest(out) %>%
  group_by(trait, significance) %>%
  mutate(p_adj = p.adjust(p = p_value, method = "bonf", n = 3),
         p_ann = case_when(p_adj <= 0.01 ~ "***",
                           p_adj <= 0.05 ~ "**",
                           p_adj <= 0.10 ~ "*",
                           # is.na(p_adj) ~ "", 
                           TRUE ~ ""),
         # Annotate the number of markers
         snp_ann = str_c("n=", nSNP)) %>%
  ungroup()
  
# Create a data.frame to plot
marker_stability_snp_prop_test_toplot <- marker_stability_snp_prop_adj %>% 
  select(trait, class, prop = null_prop) %>% 
  distinct() %>% 
  mutate(significance = "Null") %>% 
  bind_rows(marker_stability_snp_prop_adj, .) %>% 
  mutate(significance = parse_factor(significance, levels = c("Null", "Average", "Plastic", "Stable")),
         p_value = round(p_value, 4),
         p_adj = round(p_adj, 4))


# Plot
g_snp_ann_prop <- marker_stability_snp_prop_test_toplot %>% 
  ggplot(aes(x = class, y = prop, fill = significance,  group = significance)) + 
  geom_col(position = "dodge") + 
  geom_errorbar(aes(ymin = lower, ymax = upper), position = position_dodge(0.9), width = 0.5) + 
  # geom_text(aes(label = p_ann), position = position_dodge(0.9), hjust = -1, vjust = 0.5, angle = 90) +
  # geom_text(aes(label = p_adj), position = position_dodge(0.9), hjust = -1, vjust = 0.5, angle = 90) +
  # geom_text(aes(label = snp_ann), position = position_dodge(0.9), vjust = -3, size = 2) +
  scale_fill_manual(values = c("Null" = "grey75", colors), name = "Marker Type") + 
  ylab("Proportion of Markers") +
  xlab("Gene Annotation Class") +
  labs(caption = "Proximal: within 10 kbp of a gene.
       Error bars are 95% confidence interval from exact binomial test.") +
  ylim(c(0, 0.80)) +
  # facet_grid(~ trait) + 
  facet_wrap(~ trait, ncol = 2) +
  theme_poster() +
  theme(legend.position = c(0.75, 0.25),
        title = element_text(size = 10))

ggsave(filename = "marker_stability_gene_annotation.jpg", plot = g_snp_ann_prop, 
       height = 8, width = 8, path = fig_dir)












### Contribution of markers to trait variation

# Extract the unique stability coefficients
# Then add the breeding program information
S2_MET_pheno_fw_uniq <- S2_MET_pheno_mean_fw %>%
  select(trait, line_name, g, b, delta) %>% 
  distinct() %>%
  left_join(., subset(entry_list, Class == "S2TP", c(Line, Program)), 
            by = c("line_name" = "Line")) %>%
  dplyr::rename(program = Program)


## Log transform the non-linear stability estimates
S2_MET_pheno_fw_uniq_trans <- S2_MET_pheno_fw_uniq %>%
  group_by(trait) %>% 
  mutate(log_delta = log(delta)) %>%
  # Tidy
  gather(coef, value, g, b, delta, log_delta) %>% 
  filter(coef != "delta")






## Determine the proportion of variation in stability/mean that is accounted by
## all markers, then marker subsets

# Create the K matrix across all SNPs
K <- A.mat(X = M, min.MAF = 0, max.missing = 1)

# Estimate the proportion of variation due to all SNPs
mean_fw_varcomp <- S2_MET_pheno_fw_uniq_trans %>% 
  group_by(trait, coef) %>%
  do({
    # Extract data
    df <- .
    
    # Copy the line_name column
    df1 <- df %>% mutate(snp = line_name)
    
    # # Fit the model
    # fit <- relmatLmer(value ~ (1|snp), df1, relmat = list(snp = K))
    # 
    # # Extract the proportion of variance attributed to snps and to residuals
    # VarProp(fit) %>% 
    #   select(-contains("var"))
    
    fit <- mixed.solve(y = df1$value, K = K)
    data.frame(grp = c("snp", "Residual"), vcov = c(fit$Vu, fit$Ve)) %>% 
      mutate(prop = vcov / sum(vcov))
  })



# Remove extra data
mean_fw_varcomp1 <- mean_fw_varcomp %>%
  filter(grp == "snp") %>%
  select(trait, coef, prop) %>%
  ungroup()

# Plot
g_allmarker_varcomp <- mean_fw_varcomp1 %>%
  ggplot(aes(x = trait, y = prop, fill = coef)) + 
  geom_col(position = "dodge") +
  ylim(c(0,1)) + 
  theme_bw()

ggsave(filename = "allmarker_varcomp.jpg", plot = g_allmarker_varcomp, path = fig_dir,
       height = 5, width = 5, dpi = 1000)


# ## Bootstrap a confidence interval for this estimate
# mean_fw_varcomp_boot <- S2_MET_pheno_fw_uniq_trans %>%
#   group_by(trait, coef) %>%
#   do(boot_out = {
#     # Extract data
#     df <- .
# 
#     # Copy the line_name column
#     df_boots <- df %>%
#       mutate(snp = line_name) %>%
#       modelr::bootstrap(data = ., n = 100)
# 
#     # Bootstrap
#     boots_out <- df_boots$strap %>%
#       map(as.data.frame) %>%
#       map(~mixed.solve(y = .$value, K = K[.$line_name, .$line_name]))  %>%
#       map(~c(snp = .$Vu, res = .$Ve) / sum(.$Vu + .$Ve))
# 
#     # Extract the proportion of variance attributed to snps and to residuals
#     boots_out %>%
#       map(~VarProp(.) %>% select(-contains("var")))
# 
#   })
    



# Number of model fittings
n_iter <- 100

# Now for each trait, find the proportion of variation in the mean and stability
# due to different marker groups
mean_fw_martype_varcomp <- S2_MET_pheno_fw_uniq_trans %>%
  group_by(trait, coef) %>%
  do({
    
    df <- .
    
    # What is the trait we are dealing with
    tr <- unique(df$trait)
    
    # Extract the markers for the particular trait and separate by average/stable/plastic
    marker_types <- marker_mean_fw_sig %>%
      filter(trait == tr) %>% 
      split(.$significance) %>% 
      map("marker")
    
    # Use the plastic markers to make relationship matrices
    K_plas <- marker_types$Plastic %>% M[,.,drop = F] %>% A.mat(X = ., min.MAF = 0, max.missing = 1)

    # Sample size
    sample_size <- length(marker_types$Plastic)
    
    # Create random samples of the average markers
    stable_marker_samples <- rerun(.n = n_iter, sample(marker_types$Stable, size = sample_size))

    ## Fit using sommer
    y <- df$value
    X <- model.matrix(~ 1, df)
    Z_g <- model.matrix(~ -1 + line_name, df) %>%
      `colnames<-`(., colnames(K_plas))
    
    # Iterate over the samples
    var_comp_out <- stable_marker_samples %>%
      map(function(marker_sample) {
        
        # Create a relationship matrix
        K_stab <- M[,marker_sample,drop = F] %>% A.mat(X = ., min.MAF = 0, max.missing = 1)
        
        # Create a Z list
        Z <- list(
          stab = list(Z = Z_g, K = K_stab),
          plas = list(Z = Z_g, K = K_plas)
        )
        
        fit <- sommer::mmer(Y = y, X = X, Z = Z, silent = TRUE)
        
        # Extract and compute variance components
        fit$var.comp %>% 
          map_df(~unname(.))
        
      })
    
    var_comp_out %>% 
      list(., seq_along(.)) %>% 
      pmap_df(~mutate(.x, iter = .y))
    
  })


## Character vector for replacing marker types
mar_type_replace <- c("plas" = "Plastic", "stab" = "Stable")


## Summarize
# Separate variance components for stable, plastic, and average markers
mean_fw_martype_varcomp_summ <- mean_fw_martype_varcomp %>% 
  mutate(total = stab + plas  + units) %>% # Add stab + plas together
  mutate_at(vars(stab:units), funs(. / total)) %>%
  select(-units, -total) %>% 
  gather(marker_type, prop_var, stab, plas) %>%
  ungroup()

## Plot
g_mean_fw_varcomp <- mean_fw_martype_varcomp_summ %>% 
  filter(marker_type %in% names(mar_type_replace)) %>%
  mutate(marker_type = str_replace_all(marker_type, mar_type_replace)) %>%
  ggplot(aes(x = coef, y = prop_var, fill = marker_type)) + 
  geom_boxplot(position = "dodge") + 
  facet_wrap(~ trait, ncol = 2) +
  theme_bw() +
  theme(legend.position = c(0.75, 0.25))

ggsave(filename = "mean_fw_varcomp.jpg", plot = g_mean_fw_varcomp, path = fig_dir,
       height = 5, width = 5, dpi = 1000)


## Save
save_file <- file.path(result_dir, "marker_stability_varcomp.RData")
save("mean_fw_varcomp", "mean_fw_martype_varcomp", file = save_file)






## What is the relationship between marker effects for stability and marker effects
## for the mean? What about marker effect stability versus marker effects for the mean?
gs_marker_effects <- S2_MET_pheno_mean_fw %>% 
  distinct(trait, line_name, g, b) %>%
  gather(term, value, g:b) %>% 
  group_by(trait, term) %>% 
  do(data.frame(marker = colnames(M), effect = mixed.solve(y = .$value, Z = M)$u)) %>%
  spread(term, effect)

# Add the marker stability
gs_marker_effects1 <- gs_marker_effects %>% 
  left_join(., distinct(marker_mean_fw, trait, marker, a, c) %>% rename(marker_c = c, marker_a = a))

## Correlate - bootstrap for significance
marker_boot_cor <- gs_marker_effects1 %>%
  do(bootstrap(x = .$b, y = .$g, fun = "cor")) %>%
  mutate(significant = !between(0, ci_lower, ci_upper),
         cor_toprint = str_c("r = ", round(base, 3), ifelse(significant, "*", "")))



# Interesting - the correlation is stronger for GY > PH > HD
g_mean_stab_marker_eff <- gs_marker_effects1 %>%
  ggplot(aes(x = g, y = b)) +
  geom_point() +
  geom_text(data = marker_boot_cor, aes(x = Inf, y = -Inf, label = cor_toprint), vjust = -1, hjust = 1.1) + 
  facet_wrap(~trait, ncol = 2, scales = "free_x") + 
  ylab("Stability Marker Effect") +
  xlab("Genotype Mean Marker Effect") + 
  theme_bw() + theme(panel.grid = element_blank())


## Save
save_file <- file.path(result_dir, "marker_stability_varcomp.RData")
save("mean_fw_varcomp", "mean_fw_martype_varcomp", file = save_file)






## What is the relationship between marker effects for stability and marker effects
## for the mean? What about marker effect stability versus marker effects for the mean?
gs_marker_effects <- S2_MET_pheno_mean_fw %>% 
  distinct(trait, line_name, g, b) %>%
  gather(term, value, g:b) %>% 
  group_by(trait, term) %>% 
  do(data.frame(marker = colnames(M), effect = mixed.solve(y = .$value, Z = M)$u)) %>%
  spread(term, effect)

# Add the marker stability
gs_marker_effects1 <- gs_marker_effects %>% 
  left_join(., distinct(marker_mean_fw, trait, marker, a, c) %>% rename(marker_c = c, marker_a = a))

## Correlate - bootstrap for significance
marker_boot_cor <- gs_marker_effects1 %>%
  do(bootstrap(x = .$b, y = .$g, fun = "cor")) %>%
  mutate(significant = !between(0, ci_lower, ci_upper),
         cor_toprint = str_c("r = ", round(base, 3), ifelse(significant, "*", "")))



# Interesting - the correlation is stronger for GY > PH > HD
g_mean_stab_marker_eff <- gs_marker_effects1 %>%
  ggplot(aes(x = g, y = b)) +
  geom_point() +
  geom_text(data = marker_boot_cor, aes(x = Inf, y = -Inf, label = cor_toprint), vjust = -1, hjust = 1.1) + 
  facet_wrap(~trait, ncol = 2, scales = "free_x") + 
  ylab("Stability Marker Effect") +
  xlab("Genotype Mean Marker Effect") + 
  theme_bw() + theme(panel.grid = element_blank())

# Save
ggsave(filename = "mean_stab_marker_eff.jpg", plot = g_mean_stab_marker_eff,
       path = fig_dir, height = 6, width = 6, dpi = 1000)
       

## The correlation for the intercept of marker effects and marker effects for the mean
## is quote high.
## The correlation for the slope of marker effects and marker effects for the slope is
## quite high.
gs_marker_effects1 %>% 
  summarize(cor_ga = cor(g, marker_a), cor_bc = cor(b, marker_c),
            cor_gb = cor(g, b), cor_ac = cor(marker_a, marker_c))


## If the markers with the lowest and highest effect on stability are the most relevant,
## is the same true for the genotype mean?
gs_marker_effect_outliers <- gs_marker_effects %>% 
  mutate(g_upper = quantile(g, 1 - (alpha / 2)), g_lower = quantile(g, alpha / 2),
         significance = ifelse(g >= g_upper | g <= g_lower, "outlier", "normal"))
  


# Number of model fittings
n_iter <- 100

# Now for each trait, find the proportion of variation in the mean and stability
# due to different marker groups
geno_mean_martype_varcomp <- S2_MET_pheno_fw_uniq_trans %>%
  group_by(trait, coef) %>%
  do({
    
    df <- .
    
    # What is the trait we are dealing with
    tr <- unique(df$trait)
    
    # Extract the markers for the particular trait and separate by average/stable/plastic
    marker_types <- gs_marker_effect_outliers %>%
      filter(trait == tr) %>% 
      split(.$significance) %>% 
      map("marker")
    
    # Use the plastic markers to make relationship matrices
    K_out <- marker_types$outlier %>% M[,.,drop = F] %>% A.mat(X = ., min.MAF = 0, max.missing = 1)
    
    # Sample size
    sample_size <- length(marker_types$outlier)
    
    # Create random samples of the average markers
    normal_marker_samples <- rerun(.n = n_iter, sample(marker_types$normal, size = sample_size))
    
    ## Fit using sommer
    y <- df$value
    X <- model.matrix(~ 1, df)
    Z_g <- model.matrix(~ -1 + line_name, df) %>%
      `colnames<-`(., colnames(K_out))
    
    # Iterate over the samples
    var_comp_out <- normal_marker_samples %>%
      map(function(marker_sample) {
        
        # Create a relationship matrix
        K_norm <- M[,marker_sample,drop = F] %>% A.mat(X = ., min.MAF = 0, max.missing = 1)
        
        # Create a Z list
        Z <- list(
          normal = list(Z = Z_g, K = K_norm),
          outlier = list(Z = Z_g, K = K_out)
        )
        
        fit <- sommer::mmer(Y = y, X = X, Z = Z, silent = TRUE)
        
        # Extract and compute variance components
        fit$var.comp %>% 
          map_df(~unname(.))
        
      })
    
    
    # Add the iteration number and output a data.frame
    var_comp_out %>% 
      list(., seq_along(.)) %>% 
      pmap_df(~mutate(.x, iter = .y))
    
  })


## Summarize
# Separate variance components for stable, plastic, and average markers
geno_mean_martype_varcomp_summ <- geno_mean_martype_varcomp %>% 
  mutate(total = normal + outlier + units) %>% 
  mutate_at(vars(normal:units), funs(. / total)) %>%
  select(-units, -total) %>% 
  gather(marker_type, prop_var, normal, outlier) %>%
  ungroup()

## Plot
g_geno_mean_varcomp <- geno_mean_martype_varcomp_summ %>% 
  ggplot(aes(x = coef, y = prop_var, fill = marker_type)) + 
  geom_boxplot(position = "dodge") + 
  facet_wrap(~trait, ncol = 2) +
  theme_bw() +
  theme(legend.position = c(0.75, 0.25))

## plot - just the mean
mean_fw_martype_varcomp_summ %>% 
  filter(coef == "g") %>%
  ggplot(aes(x = trait, y = prop_var, fill = marker_type)) + 
  geom_boxplot(position = "dodge") + 
  theme_bw()

ggsave(filename = "geno_mean_varcomp.jpg", plot = g_geno_mean_varcomp, path = fig_dir,
       height = 5, width = 5, dpi = 1000)






##### Marker stability as effect estimates for phenotypic stability
# Use the marker stability estimates as marker effects - predict phenotypic stability

# Create matrices with these "effects"
marker_stability_effect <- marker_mean_fw %>% 
  distinct(trait, marker, b) %>% 
  split(.$trait) %>% 
  map(~select(., -trait) %>% as.data.frame() %>% column_to_rownames("marker") %>% as.matrix())

# Split the phenotypic stability estimates
pheno_mean_fw_split <- S2_MET_pheno_mean_fw %>% 
  distinct(trait, line_name, b) %>% 
  split(.$trait)

## Predictions of stability
marker_stability_pred_acc <- list(marker_stability_effect, pheno_mean_fw_split) %>%
  pmap(~mutate(.y, marker_b_pred = (M %*% (.x)) + 1)) %>%
  bind_rows()

# Accuracy
marker_stability_pred_acc %>% 
  group_by(trait) %>% 
  summarize(marker_b_pred_acc = cor(b, marker_b_pred))

### Pretty high!



## Predictions of stability
marker_stability_pred_acc <- list(marker_stability_effect, pheno_mean_fw_split) %>%
  pmap(~mutate(.y, marker_b_pred = (M %*% (.x)) + 1)) %>%
  bind_rows()

# Accuracy
marker_stability_pred_acc %>% 
  group_by(trait) %>% 
  summarize(marker_b_pred_acc = cor(b, marker_b_pred))

### Pretty high!

# Plot
g_marker_stability_pred <- marker_stability_pred_acc %>% 
  ggplot(aes(x = b, y = marker_b_pred)) +
  geom_point() +
  facet_wrap(~ trait, ncol = 2) + 
  xlab("Linear Phenotypic Stability") +
  ylab("Predicted Linear Stability Using Marker Stability") +
  theme_bw()

ggsave(filename = "marker_stability_pred.jpg", plot = g_marker_stability_pred, path = fig_dir,
       height = 6, width = 8, dpi = 1000)


# Calculate marker effects for estimates and stability and compare them to the
# estimates of marker effect stability
pheno_stab_marker_effect <- pheno_mean_fw_split %>%
  map(~data.frame(trait = .$trait[1], pheno_stab_marker_effect = mixed.solve(y = .$b, Z = M)$u))

## Stability marker effects (sme) versus marker effect stability (mes)
sme_vs_mes <- list(marker_stability_effect, pheno_stab_marker_effect) %>%
  pmap_df(~cbind(.y, .x))


# Calculate marker effects for estimates and stability and compare them to the
# estimates of marker effect stability
pheno_stab_marker_effect <- pheno_mean_fw_split %>%
  map(~data.frame(trait = .$trait[1], pheno_stab_marker_effect = mixed.solve(y = .$b, Z = M)$u))

## Stability marker effects (sme) versus marker effect stability (mes)
sme_vs_mes <- list(marker_stability_effect, pheno_stab_marker_effect) %>%
  pmap_df(~cbind(.y, .x))

# Correlation
sme_vs_mes %>% 
  group_by(trait) %>%
  summarize(cor = cor(pheno_stab_marker_effect, b))

# Plot
g_sme_vs_mes <- sme_vs_mes %>% 
  ggplot(aes(x = pheno_stab_marker_effect, y = b)) +
  geom_point() +
  facet_wrap(~ trait, ncol = 2) + 
  xlab("Phenotypic Linear Stability Marker Effect") +
  ylab("Marker Effect Linear Stability") +
  theme_bw()



# Plot
g_sme_vs_mes <- sme_vs_mes %>% 
  ggplot(aes(x = pheno_stab_marker_effect, y = b)) +
  geom_point() +
  facet_wrap(~ trait, ncol = 2) + 
  xlab("Phenotypic Linear Stability Marker Effect") +
  ylab("Marker Effect Linear Stability") +
  theme_bw()

ggsave(filename = "marker_effect_stability_linear_stability.jpg", plot = g_sme_vs_mes, path = fig_dir,
       height = 6, width = 8, dpi = 1000)



### If marker effects stability acts additively on phenotypic stability, then
### we could look at different chromosomes to construct more stable genotypes
# Get the allele stability states for all genotypes
marker_stability_states <- marker_stability_effect %>% 
  map(~M * c(. - 1)) %>%
  map(~as.data.frame(.) %>% rownames_to_column("line_name") %>% gather(marker, effect_state, -line_name)) %>%
  list(., names(.)) %>%
  pmap_df(~mutate(.x, trait = .y))

# Combine with marker metadata
marker_stability_states1 <- marker_stability_states %>% 
  left_join(., snp_info) %>%
  as_data_frame()

# For each line, sum the effects of each chromosome
marker_stability_states_sum <- marker_stability_states1 %>%
  group_by(trait, line_name, chrom) %>% 
  summarize(chrom_sum = sum(effect_state))

marker_stability_states_sum %>%
  filter(trait == trait[1]) %>%
  # Sort on line breeding program
  left_join(., select(entry_list, line_name = Line, program = Program)) %>%
  ggplot(aes(x = line_name, y = chrom, fill = chrom_sum)) +
  geom_tile() +
  scale_fill_gradient2(low = "red", mid = "white", high = "blue") +
  facet_grid(~ program, scales = "free_x", space = "free_x") +
  theme_bw() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        panel.spacing.x = unit(0, "lines"))











# ######
# # Appendix
# ######
# 
# 
# 
# #### Overlap with Stability and Mean Loci
# ## Also look at previously described QTL
# 
# # Load the GWAS results for genotypic mean and phenotypic stability and see if the 
# # markers with the highest plasticity/stability measure overlap
# 
# load(file.path(result_dir, "gwas_adjusted_significant_results.RData"))
# 
# # Read in the QTL metadata files that were downloaded from T3
# qtl_meta <- map_df(list.files(data_dir, pattern = "meta", full.names = TRUE), read_csv)
# # Read in association data from Pauli2014 and Wang2012
# bcap_association <- read_csv(file.path(data_dir, "BCAP_association_qtl.csv")) %>%
#   mutate(position = parse_number(position)) %>%
#   select(trait, marker, chrom = chromosome, pos = position, gene:feature, reference) %>%
#   filter(trait %in% c("GrainYield", "HeadingDate", "PlantHeight"))
# 
# 
# # Rename the traits
# qtl_meta_df <- qtl_meta %>%
#   mutate(trait = case_when(trait == "grain yield" ~ "GrainYield",
#                            trait == "heading date" ~ "HeadingDate",
#                            trait == "plant height" ~ "PlantHeight"))
# 
# # Remove unmapped QTL and change the chromosome name
# qtl_meta_use <- qtl_meta_df %>% 
#   filter(!chromosome %in% c("chrUNK", "chrUn")) %>% 
#   mutate(chrom = as.integer(parse_number(chromosome))) %>%
#   select(trait, marker, chrom, pos = position, gene, feature) %>%
#   arrange(trait, chrom, pos)
# 
# # Combine data
# qtl_meta_use1 <- bind_rows(qtl_meta_use, bcap_association)
# 
# 
# 
# ### First highlight some examples for known genes / QTL for heading date and plant height
# 
# ## Heading Date
# tr <- "HeadingDate"
# chr <- 5
# start <- 580e6
# end <- 620e6
# 
# 
# ## Add information on VRN
# vrn1_data <- data.frame(chrom = 5, gene_start = 599135017, gene_end = 599147377, gene_id = "Vrn-H1")
# 
# # Get the QTL information
# hd_qtl <- qtl_meta_use1 %>% 
#   filter(trait == tr, chrom == chr)
# 
# 
# 
# ## Subset the chromosome of interest
# g_hd_marstab <- S2_MET_marker_eff_pheno_fw_sig %>%   
#   filter(coef == "b", trait == tr, chrom == chr, between(pos, start, end)) %>%
#   ggplot(aes(x = pos / 1000000, y = estimate, col = significance)) + 
#   geom_point() + 
#   geom_hline(aes(yintercept = lower_perc), lty = 2) +
#   geom_hline(aes(yintercept = upper_perc), lty = 2) +
#   scale_color_manual(values = colors, name = NULL) +
#   ylab("Marker Effect\nPlasticity Estimate") +
#   xlab("Position (Mbp)") +
#   labs(title = "Heading Date") + 
#   theme_poster() +
#   theme(legend.position = c(0.15, 0.1),
#         axis.title.x = element_blank(),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         panel.border = element_blank())
# 
# 
# ## Create a separate gene annotation plot
# g_hd_ann <- vrn1_data %>%
#   ggplot(aes(y = 1, yend = 1)) +
#   xlim(c(start / 1e6, end / 1e6)) + 
#   geom_segment(aes(x = (gene_start - 100000) / 1000000, xend = (gene_end + 100000) / 1000000), lwd = 5) +
#   geom_text(aes(x = gene_start / 1000000, label = gene_id), vjust = 2, fontface = "italic") +
#   # Add significant QTL
#   geom_segment(data = hd_qtl, aes(x = (pos - 100000) / 1000000, xend = (pos + 100000) / 1000000), lwd = 5) +
#   xlab("Position on Chromosome 5 (Mbp)") +
#   ylab("Known\nGenes/QTL") +
#   ylim(c(0, 1)) +
#   theme_poster() +
#   theme(axis.line.x = element_line(),
#         axis.text.y = element_blank(),
#         # axis.title.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         panel.border = element_blank(),
#         panel.grid = element_blank())
# 
# ## Add the GWAS analysis of stabililty
# g_hd_stab_gwas <- gwas_pheno_mean_fw_adj %>% 
#   filter(trait == tr, coef == "b", chrom == chr, between(pos, start, end)) %>%
#   ggplot(aes(x = pos / 1000000, y = qvalue_neg_log10)) + 
#   geom_point() +
#   ylab(expression(atop(-log[10](italic(q)),'Phenotypic Plasticty'))) +
#   theme_poster() +
#   theme(axis.title.x = element_blank(),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         panel.border = element_blank())
# 
# 
# ## Combine
# g_hd_example <- plot_grid(
#   g_hd_marstab, 
#   g_hd_stab_gwas,
#   g_hd_ann,
#   ncol = 1, rel_heights = c(0.8, 0.4, 0.25), align = "hv", axis = "lr")
# 
# ## Save
# ggsave(filename = "hd_marstab_annotation_poster.jpg", plot = g_hd_example, path = fig_dir,
#        height = 8, width = 6, dpi = 1000)
# 
# 
# ## Plant Height
# tr <- "PlantHeight"
# chr <- 6
# start <- 0
# end <- 20e6
# 
# 
# ph_qtl <- qtl_meta_use1 %>% 
#   filter(trait == tr, chrom == chr)
# 
# ## Plot of marker effect plasticity
# g_ph_marstab <- S2_MET_marker_eff_pheno_fw_sig %>%   
#   filter(coef == "b", trait == tr, chrom == chr, between(pos, start, end)) %>%
#   ggplot(aes(x = pos / 1000000, y = estimate, col = significance)) + 
#   geom_point() + 
#   geom_hline(aes(yintercept = lower_perc), lty = 2) +
#   geom_hline(aes(yintercept = upper_perc), lty = 2) +
#   scale_color_manual(values = colors, name = NULL) +
#   ylab("Marker Plasticity Estimate") +
#   xlab("Position (Mbp)") +
#   theme_poster() +
#   theme(legend.position = c(0.15, 0.05),
#         axis.title.x = element_blank(),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         panel.border = element_blank())
# 
# 
# ## Create a separate gene annotation plot
# g_ph_ann <- ph_qtl %>%
#   ggplot(aes(y = 0, yend = 0)) +
#   xlim(c(start / 1e6, end / 1e6)) + 
#   geom_segment(aes(x = (pos - 100000) / 1000000, xend = (pos + 100000) / 1000000), lwd = 5) +
#   xlab("Position on Chromosome 6 (Mbp)") +
#   theme_poster() +
#   theme(axis.line.x = element_line(),
#         axis.text.y = element_blank(),
#         axis.title.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         panel.border = element_blank(),
#         panel.grid = element_blank())
# 
# ## Add the GWAS analysis of stabililty
# g_ph_stab_gwas <- gwas_pheno_mean_fw_adj %>% 
#   filter(trait == tr, coef == "b", chrom == chr, between(pos, start, end)) %>%
#   ggplot(aes(x = pos / 1000000, y = qvalue_neg_log10)) + 
#   geom_point() +
#   ylab(expression(-log[10](italic(q))~'Phenotypic Plasticty')) +
#   theme_poster() +
#   theme(axis.title.x = element_blank(),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         panel.border = element_blank())
# 
# ## Combine
# g_ph_example <- plot_grid(
#   g_ph_marstab, 
#   g_ph_stab_gwas,
#   g_ph_ann,
#   ncol = 1, rel_heights = c(0.8, 0.4, 0.2), align = "hv", axis = "lr")
# 
# ## Save
# ggsave(filename = "ph_marstab_annotation_poster.jpg", plot = g_ph_example, path = fig_dir,
#        height = 8, width = 5, dpi = 1000)
# 
# 
# ## Combine HD and PH
# g_hd_ph_example <- plot_grid(g_hd_example, g_ph_example, ncol = 2)
# ## Save
# ggsave(filename = "hd_ph_marstab_annotation_poster.jpg", plot = g_hd_ph_example, path = fig_dir,
#        height = 8, width = 10, dpi = 1000)
# 
# ## Plant Height
# tr <- "PlantHeight"
# chr <- 6
# start <- 0
# end <- 20e6
# 
# 
# ph_qtl <- qtl_meta_use1 %>% 
#   filter(trait == tr, chrom == chr)
# 
# ## Plot of marker effect plasticity
# g_ph_marstab <- S2_MET_marker_eff_pheno_fw_sig %>%   
#   filter(coef == "b", trait == tr, chrom == chr, between(pos, start, end)) %>%
#   ggplot(aes(x = pos / 1000000, y = estimate, col = significance)) + 
#   geom_point() + 
#   geom_hline(aes(yintercept = lower_perc), lty = 2) +
#   geom_hline(aes(yintercept = upper_perc), lty = 2) +
#   scale_color_manual(values = colors, name = NULL) +
#   ylab("Marker Plasticity Estimate") +
#   xlab("Position (Mbp)") +
#   theme_poster() +
#   theme(legend.position = c(0.15, 0.05),
#         axis.title.x = element_blank(),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         panel.border = element_blank())
# 
# 
# ## Create a separate gene annotation plot
# g_ph_ann <- ph_qtl %>%
#   ggplot(aes(y = 0, yend = 0)) +
#   xlim(c(start / 1e6, end / 1e6)) + 
#   geom_segment(aes(x = (pos - 100000) / 1000000, xend = (pos + 100000) / 1000000), lwd = 5) +
#   xlab("Position on Chromosome 6 (Mbp)") +
#   theme_poster() +
#   theme(axis.line.x = element_line(),
#         axis.text.y = element_blank(),
#         axis.title.y = element_blank(),
#         axis.ticks.y = element_blank(),
#         panel.border = element_blank(),
#         panel.grid = element_blank())
# 
# ## Add the GWAS analysis of stabililty
# g_ph_stab_gwas <- gwas_pheno_mean_fw_adj %>% 
#   filter(trait == tr, coef == "b", chrom == chr, between(pos, start, end)) %>%
#   ggplot(aes(x = pos / 1000000, y = qvalue_neg_log10)) + 
#   geom_point() +
#   ylab(expression(-log[10](italic(q))~'Phenotypic Plasticty')) +
#   theme_poster() +
#   theme(axis.title.x = element_blank(),
#         axis.text.x = element_blank(),
#         axis.ticks.x = element_blank(),
#         panel.border = element_blank())
# 
# ## Combine
# g_ph_example <- plot_grid(
#   g_ph_marstab, 
#   g_ph_stab_gwas,
#   g_ph_ann,
#   ncol = 1, rel_heights = c(0.8, 0.4, 0.2), align = "hv", axis = "lr")
# 
# ## Save
# ggsave(filename = "ph_marstab_annotation_poster.jpg", plot = g_ph_example, path = fig_dir,
#        height = 8, width = 5, dpi = 1000)
# 
# 
# ## Combine HD and PH
# g_hd_ph_example <- plot_grid(g_hd_example, g_ph_example, ncol = 2)
# ## Save
# ggsave(filename = "hd_ph_marstab_annotation_poster.jpg", plot = g_hd_ph_example, path = fig_dir,
#        height = 8, width = 10, dpi = 1000)
# 
# 
# 
# 
# ## Plot examples of the stability of genotypes/markers
# pheno_stable_ex <- S2_MET_pheno_mean_fw %>% 
#   group_by(trait) %>% 
#   filter(b == min(b) | b == max(b) | abs(b - 1) == min(abs(b - 1))) %>%
#   mutate(class = case_when(
#     b == max(b) ~ "Plastic",
#     b == min(b) ~ "Stable",
#     TRUE ~ "Average")) %>%
#   ungroup()
# 
# g_pheno_stable_ex <- pheno_stable_ex %>% 
#   # mutate(value = value - h) %>%
#   ggplot(aes(x = h, y = value, color = class)) + 
#   geom_point() + 
#   geom_smooth(method = "lm", se = FALSE) +  
#   scale_color_manual(values = colors, name = "Response Type") +
#   facet_wrap(~ trait, ncol = 1, scales = "free") +
#   ylab(expression("Phenotypic Value "~(italic(y[ij])))) + 
#   xlab(expression("Environmental Effect "~(italic(t[j])))) +
#   labs(title = "Phenotypic Plasticity") + 
#   theme_poster() +
#   theme(legend.position = "bottom", 
#         legend.text = element_text(size = 14),
#         title = element_text(size = 16))
# 
# 
# mark_stable_ex <- S2_MET_marker_mean_fw %>%
#   group_by(trait) %>% 
#   filter(b == min(b) | b == max(b) | abs(b - 0) == min(abs(b - 0))) %>%
#   mutate(class = case_when(
#     b == max(b) ~ "Plastic",
#     b == min(b) ~ "Stable",
#     TRUE ~ "Average")) %>%
#   ungroup()
# 
# g_mark_stable_ex <- mark_stable_ex %>% 
#   # mutate(effect = effect + h) %>%
#   ggplot(aes(x = h, y = effect, color = class)) + 
#   geom_point() + 
#   geom_smooth(method = "lm", se = FALSE) +  
#   scale_color_manual(values = colors, name = NULL) +
#   facet_wrap(~ trait, ncol = 1, scales = "free") +
#   ylab(expression("Marker Effect "~(italic(alpha[jp])))) + 
#   xlab(expression("Environmental Effect "~(italic(t[j])))) +
#   labs(title = "Marker Effect Plasticity",
#        caption = expression("Environmental effect"~(italic(t[j]))~"not added to marker effect to highlight reponse.")) + 
#   theme_poster() +
#   theme(legend.position = "none",
#         plot.title = element_text(size = 16),
#         plot.caption = element_text(size = 10))
# 
# ## Combine
# g_stable_ex <- plot_grid(g_pheno_stable_ex + theme(legend.position = "none"), 
#                          g_mark_stable_ex, ncol = 2, align = "hv")
# g_stable_ex1 <- plot_grid(g_stable_ex, get_legend(g_pheno_stable_ex), ncol = 1, rel_heights = c(0.95, 0.05))
# 
# ## Save
# ggsave(filename = "stability_example_poster.jpg", plot = g_stable_ex1, path = fig_dir,
#        height = 10, width = 8, dpi = 1000)
# 
# 
# 
