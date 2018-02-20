## S2MET Mapping
## 
## Marker Stability Analysis
## 


library(tidyverse)
library(stringr)
library(readxl) 
library(GenomicRanges)
library(qvalue)
library(lme4qtl)
library(patchwork)
library(neyhart)
library(rrBLUP)

repo_dir <- getwd()

# Project and other directories
source("C:/Users/Jeff/Google Drive/Barley Lab/Projects/S2MET_Mapping/source.R")

# Load the stability results
load(file.path(result_dir, "S2MET_pheno_mean_fw_results.RData"))

# # List the files containing marker effect by environment estimates
# mxe_files <- list.files(result_dir, pattern = "marker_by_env", full.names = TRUE)
# 
# # Iterate over files and save in a list
# mxe_list <- list()
# for (file in mxe_files) {
#   load(file)
#   trait <- str_extract(string = file, "[A-Za-z]*.RData") %>% str_replace(pattern = ".RData", "")
#   
#   mxe_list[[trait]] <- marker_score_out
# }

# Load the mxe results
load(file.path(result_dir, "S2MET_marker_ranef_by_env.RData"))

# Get the marker information
snp_info <- S2TP_imputed_multi_genos_hmp %>%
  select(marker = rs, chrom, pos, cM_pos) %>%
  # Correct the cM position for BOPA snps
  mutate(cM_pos = if_else(str_detect(marker, "^S"), cM_pos, cM_pos / 1000))

# Filter the BLUEs to use
S2_MET_BLUEs_use <- S2_MET_BLUEs %>% 
  filter(line_name %in% c(tp_geno)) %>%
  droplevels()

M <- S2TP_imputed_multi_genos_mat[tp_geno,]




## Manage the mxe data and convert to a tidy data.frame

# # Tidy up
# marker_by_environment <- mxe_list %>%
#   map(~map_df(., ~list(., names(.)) %>% pmap_df(~mutate(.x, marker = .y)))) %>%
#   list(., names(.)) %>%
#   pmap_df(~mutate(.x, trait = .y)) %>%
#   as_data_frame() %>%
#   select(trait, marker, names(.))
# 
# # Rename
# mxe_df <- marker_by_environment %>% 
#   mutate(environment = str_extract(term, "[A-Z]{3}[0-9]{2}")) %>%
#   select(trait, marker, environment, estimate)

mxe_df <- marker_by_env %>% 
  unnest(marker_effect) %>% 
  ungroup()

# Combine with the environmental mean
mxe_df1 <- mxe_df %>%
  left_join(., distinct(S2_MET_pheno_mean_fw, trait, environment, h), 
            by = c("environment", "trait"))

# Calculate the marker effect stability
marker_effect_stab <- mxe_df1 %>%
  group_by(trait, marker) %>%
  do({
    df <- . 
    # Fit the model
    fit <- lm(effect ~ h, df)
    data_frame(b = coef(fit)[2], 
               b_std_error = summary(fit)$coefficients[2,"Std. Error"],
               delta = mean(resid(fit)^2),
               df = df.residual(fit),
               fit = list(fit)) })
    


# Add the snp information
S2_MET_marker_mean_fw <- mxe_df1 %>% 
  full_join(., marker_effect_stab, by = c("trait", "marker")) %>% 
  left_join(., snp_info, by = "marker") %>% 
  select(marker, chrom:cM_pos, trait, environment, effect, h:delta, df)

# Save the data
save_file <- file.path(result_dir, "S2MET_marker_mean_fw_results.RData")
save("S2_MET_marker_mean_fw", "marker_effect_stab", file = save_file)


## Plot the distribution of stability estimates

## Transform
S2_MET_marker_mean_fw_trans <- S2_MET_marker_mean_fw %>% 
  distinct(marker, chrom, pos, trait, b, delta) %>%
  mutate(log_delta = log(delta), b = b + 1)


## Plot the stability terms individually, then combine
# First define a common list of ggplot modifiers
g_mod <- list(
  geom_density(aes(fill = "blue")),
  xlab("Estimate"),
  theme_bw() +
  theme(axis.title.y = element_blank()) )

# Just plot linear stability
g_marker_fw_dens_b <- S2_MET_marker_mean_fw_trans %>%
  ggplot(aes(x = b)) + 
  labs(title = "Linear Stability") +
  facet_wrap( ~ trait, ncol = 1) +
  scale_fill_brewer(palette = "Set2", guide = FALSE) +
  g_mod


# Just plot non-linear stability
g_marker_fw_dens_delta <- S2_MET_marker_mean_fw_trans %>%
  ggplot(aes(x =  log_delta)) + 
  labs(title = "Non-Linear Stability") +
  scale_fill_brewer(palette = "Set2", guide = FALSE) +
  facet_wrap( ~ trait, ncol = 1, scales = "free_x") +
  g_mod

# Add the plots together
g_marker_fw_dist <- g_marker_fw_dens_b + g_marker_fw_dens_delta

save_file <- file.path(fig_dir, "marker_stability_estimate_distriubtions.jpg")
ggsave(filename = save_file, plot = g_marker_fw_dist, height = 5.6, width = 5, dpi = 500)




# First group markers into those that are highly stable, highly sensitive, or neither using empirical thresholds

# Tidy up
S2_MET_marker_mean_fw_tidy <- S2_MET_marker_mean_fw_trans %>% 
  select(-delta) %>% 
  gather(coef, estimate, -marker:-trait)
  

# What should be the cutoff level?
alpha <- 0.05

# For each trait, calculate empirical thresholds for significance
S2_MET_marker_eff_pheno_fw_sig <- S2_MET_marker_mean_fw_tidy %>%
  filter(coef == "b") %>% 
  group_by(trait) %>% 
  # mutate(estimate = scale(estimate)) %>%
  mutate(lower_perc = quantile(estimate, alpha / 2), 
         upper_perc = quantile(estimate, 1 - (alpha / 2))) %>%
  ungroup() %>%
  mutate(significance = case_when(estimate >= upper_perc ~ "Plastic",
                                  estimate <= lower_perc ~ "Stable",
                                  TRUE ~ "Average"),
         marker_type = if_else(str_detect(marker, "^S"), "GBS", "BOPA"))

# Calculate a table of GBS markers and sensitive/stable/notsignificant
mar_stab_table <- S2_MET_marker_eff_pheno_fw_sig %>% 
  group_by(trait, marker_type, significance) %>% 
  summarize(n = n()) %>% 
  mutate(prop = n / sum(n))

mar_stab_table %>% 
  select(-n) %>% 
  spread(significance, prop)


# Plot the distributions of marker effect stabilities and show the empirical threholds
g_b_dist <- S2_MET_marker_eff_pheno_fw_sig %>% 
  ggplot(aes(x = estimate, fill = "red")) + 
  geom_density() + 
  geom_vline(aes(xintercept = upper_perc)) +
  geom_vline(aes(xintercept = lower_perc)) + 
  facet_wrap(~ trait, ncol = 1, scales = "free_x") + 
  scale_fill_discrete(guide = FALSE) + 
  theme_bw()

# Boxplots
S2_MET_marker_mean_fw_tidy %>% 
  mutate(coef = if_else(coef == "b", "linear", "non-linear")) %>% 
  ggplot(aes(x = trait, y = estimate)) + 
  geom_boxplot() +
  facet_wrap(~ coef, scales = "free_y") +
  theme_bw()

# Calculate the mean stability per group
S2_MET_marker_eff_pheno_fw_sig %>% 
  group_by(trait, significance) %>% 
  summarize(mean = mean(estimate)) %>% 
  spread(significance, mean)



# Plot the stability estimates across the genome with the empirical threshold

# Color theme
colors <- c("black", umn_palette(n = 4)[-c(1:2)]) %>%
  set_names(c("Average", "Plastic", "Stable"))


g_mar_stab <- S2_MET_marker_eff_pheno_fw_sig %>%   
  ggplot(aes(x = pos / 1000000, y = estimate, col = significance)) + 
  geom_point() + 
  geom_hline(aes(yintercept = lower_perc), lty = 2) +
  geom_hline(aes(yintercept = upper_perc), lty = 2) +
  scale_color_manual(values = colors, name = "Marker Group") +
  ylab("Linear Stability Estimate") +
  xlab("Position (Mbp)") +
  facet_grid(trait ~ chrom, scales = "free", space = "free_x", switch = "x") +
  theme_manhattan()


# Save
save_file <- file.path(fig_dir, "marker_stability_genome.jpg")
ggsave(filename = save_file, plot = g_mar_stab, height = 6.5, width = 9, dpi = 500)



## Determine the proportion of GxE variance ($V_{GE}$) that is explained by 
## sensitive/stable markers vs not signficant markers

# Load the results
load(file.path(result_dir, "S2MET_pheno_mar_fw_varcomp_GrainYield.RData"))

var_comp_prop <- var_comp_out$GrainYield %>% 
  bind_rows() %>%
  by_row(sum, .collate = "cols", .to = "tot_var") %>%
  mutate_at(vars(-tot_var), funs(prop = . / tot_var))
  
# Test if the observations come from the same distribution
var_comp_prop %>% 
  summarize(ks_test = ks.test(x = gt_s_prop, gt_ns_prop)$p.value)
  
# Plot density plots for each type of marker
var_comp_prop %>%
  select(gt_s_prop, gt_ns_prop) %>%
  gather(term, prop) %>%
  ggplot(aes(x = prop, fill = term)) +
  geom_density(alpha = 0.75) +
  theme_bw()





## Gene Annotation

# Look at gene annotation with relation to the different groups of SNPs. 
# Determine the proportion of SNPs in genes, gene-proximal, or intergenic


## Load the barley annotation
ann_dir <- "C:/Users/Jeff/Google Drive/Barley Lab/Projects/Genomics/Annotation"
load(file.path(ann_dir, "barley_genomic_ranges.RData"))

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
  full_join(., select(S2_MET_marker_eff_pheno_fw_sig, marker, trait, significance),
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

## Perform an binomial exact test
## H0: the proportion of SNPs in each significance group in each class is
## not different than the proportion across all SNPs
marker_stability_snp_prop_test <- marker_stability_snp_prop %>% 
  group_by(trait, significance, class) %>%
  mutate(binom_test = list({binom.test(x = nSNP, n = significance_SNP, p = null_prop)})) %>%
  ungroup() %>% 
  mutate(p_value = map_dbl(binom_test, "p.value"),
         # Adjust p-values to reflect multiple testing (3 tests per marker group)
         p_adj =  map_dbl(p_value, ~p.adjust(p = ., method = "bonf", n = 3)),
         p_ann = case_when(p_adj <= 0.01 ~ "***",
                           p_adj <= 0.05 ~ "**",
                           p_adj <= 0.10 ~ "*",
                           # is.na(p_adj) ~ "",
                           TRUE ~ ""))
  
# Create a data.frame to plot
marker_stability_snp_prop_test_toplot <- marker_stability_snp_prop_test %>% 
  select(trait, class, prop = null_prop) %>% 
  distinct() %>% 
  mutate(significance = "Null") %>% 
  bind_rows(marker_stability_snp_prop_test, .) %>% 
  mutate(significance = parse_factor(significance, levels = c("Null", "Average", "Plastic", "Stable")),
         p_value = round(p_value, 4),
         p_adj = round(p_adj, 4))


# Plot
g_snp_ann_prop <- marker_stability_snp_prop_test_toplot %>% 
  ggplot(aes(x = class, y = prop, fill = significance,  group = significance)) + 
  geom_col(position = "dodge") + 
  geom_text(aes(label = p_ann), position = position_dodge(0.9), hjust = -1, vjust = 0.5, angle = 90) +
  # geom_text(aes(label = p_adj), position = position_dodge(0.9), hjust = -1, vjust = 0.5, angle = 90) + 
  scale_fill_manual(values = rev(umn_palette(n = 6)[-c(1:2)]), name = "Marker Group") + # Add custom colors
  ylab("Proportion of Markers") +
  xlab("Annotation Class") +
  ylim(c(0, 1)) +
  facet_grid(~ trait) + 
  theme_bw()

save_file <- file.path(fig_dir, "marker_stability_gene_annotation.jpg")
ggsave(filename = save_file, plot = g_snp_ann_prop, height = 5, width = 8)




#### Overlap with Stability and Mean Loci

# Load the GWAS results for genotypic mean and phenotypic stability and see if the 
# markers with the highest plasticity/stability measure overlap

load(file.path(result_dir, "S2MET_gwas_sig_loci.RData"))

# First determine the markers with the highest stability/plasticity per trait per chromosome
marker_stab_outliers <- S2_MET_marker_eff_pheno_fw_sig %>% 
  filter(significance != "Average")
  # group_by(trait, coef, chrom) %>% 
  # filter(estimate == max(estimate) | estimate == min(estimate))

## Add MAF and program specific allele frequency information
# Calculate the allele frequency of the 1 allele
af1 <- {colMeans(M + 1) / 2} %>%
  data.frame(marker = names(.), af = ., row.names = NULL, stringsAsFactors = FALSE)

# Subset the marker matrix for these SNPs and count the number of lines
# from each program with each allele
marker_stab_outliers_prop <- marker_stab_outliers %>%
  left_join(., af1, by = "marker") %>% 
  group_by(trait, marker) %>%
  do({
    df <- .
    
    # Extract the data on the marker
    marker_geno <- M[,df$marker, drop = FALSE] %>% 
      as.data.frame() %>% 
      rownames_to_column("line_name") %>% 
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

# Create a table to print
marker_stab_outliers_toprint <- marker_stab_outliers_prop %>%
  ungroup() %>%
  left_join(., snp_info) %>%
  select(trait, coef, marker, chrom, pos, cM = cM_pos, estimate, significance, af, AB:WA) %>%
  mutate(coef = str_replace_all(coef, c("b" = "Linear Stability", "log_delta" = "Non-Linear Stability", 
                                  "g" = "Genotype Mean")),
         MAF = pmin(af, 1 - af)) %>%
  rename_at(.vars = vars(trait:pos), .funs = str_to_title) %>%
  select(Trait,Type = significance, Marker:cM, MAF, Estimate = estimate, AB:WA) %>%
  arrange(Trait, Type, Chrom, Pos)

# Write the table
save_file <- file.path(fig_dir, "marker_stability_outliers.csv")
write_csv(x = marker_stab_outliers_toprint, path = save_file)



## Further analysis


# Create Grange objects
marker_stab_outliers_grange <- marker_stab_outliers %>% 
  group_by(trait, significance) %>% 
  do(queryGrange = makeGRangesFromDataFrame(df = ., keep.extra.columns = TRUE, start.field = "pos", end.field = "pos"))

# Create granges with the signficant GWAS SNPs
sig_markers_grange <- sig_markers_grange_df %>% 
  group_by(trait, coef) %>% 
  do(subjectGrange = GRanges(.))

# Intersect with the multi-locus mixed model results
grange_combined <- inner_join(x = sig_markers_grange, y = marker_stab_outliers_grange) %>%
  ungroup()


## Group by row, and find any overlap between the stable/plastic markers and QTL for the
## specific trait/coef
grange_combined_overlap <- grange_combined %>%
  rowwise() %>%
  do({
    df <- .
    
    # Calculate the distance from each of the markers to the significant loci (on the 
    # appropriate chromosome )
    hits <- distanceToNearest(x = df$queryGrange, subject = df$subjectGrange)
    
    # Filter for overlap and get the granges for the subject and query
    hits_grange <- hits %>% 
      as_data_frame() %>% 
      filter(distance == 0) %>%
      mutate(queryGrange = map(queryHits, ~df$queryGrange[.]),
             subjectGrange = map(subjectHits, ~df$subjectGrange[.]))
    
    # Combine with the original df
    hits_grange %>% 
      mutate(trait = df$trait, coef = df$coef, significance = df$significance) %>% 
      select(trait:significance, distance, queryGrange, subjectGrange)
    
  })    

# Convert the Granges to df
grange_combined_overlap_df <- grange_combined_overlap %>% 
  ungroup() %>% 
  mutate_at(vars(contains("Grange")), funs(map(., as.data.frame)))



# Find the proprotion of markers that overlapped for each trait
# 232 is the number of plastic or stable markers and 2 is the number of classes (plastic and stable)
# This tells us the proportion of stable/plastic markers that overlapped with any
# of the significant regions for
grange_combined_overlap_df %>% 
  mutate(query_markers = map_chr(queryGrange, "marker")) %>% 
  group_by(trait, coef) %>% 
  summarize(prop_query_markers = n_distinct(query_markers) / (232 * 2))



# Convert to a more readable DF
grange_combined_overlap_df1 <- grange_combined_overlap_df %>% 
  unnest() %>% 
  select(trait, coef, significance, chrom = seqnames, query_pos = start, 
         query_marker = marker, query_estimate = estimate, 
         query_marker_type = marker_type, subject_start = start1, subject_end = end1, 
         subject_marker = marker1, alpha, subject_pos = snp_pos)

## Summarize per trait and signficant associations
grange_combined_overlap_summary <- grange_combined_overlap_df1 %>% 
  group_by(trait, coef, subject_marker, significance) %>% 
  summarize(chrom = unique(chrom), 
            n_marker = n(), 
            any_perf_marker = any(query_marker == subject_marker),
            avg_distance = mean(abs(query_pos - subject_pos)))



# Calculate the distance and LD between the SNPs
grange_combined_overlap_summ <- grange_combined_overlap_df1 %>% 
  rowwise() %>%
  do({
    df <- .
    
    distance <- abs(df$query_pos - df$subject_pos)
    LD <- (cor(M[,c(df$query_marker, df$subject_marker)])^2)[1,2]
    
    df %>% 
      as.data.frame() %>%
      select(trait:chrom) %>% 
      mutate(distance = distance, LD = LD)
    
  })





## Contribution of markers to trait variation


# Create the K matrix across all SNPs
K <- A.mat(X = S2TP_imputed_multi_genos_mat, min.MAF = 0, max.missing = 1)

## Estimate heritability by using the K matrix to model genotype effect (i.e. SNP heritability)
S2_MET_fw_herit <- S2_MET_pheno_fw_uniq_trans %>% 
  group_by(trait, term) %>%
  do({
    # Extract data
    df <- .
    
    # Copy the line_name column
    df1 <- df %>% mutate(snp = line_name)
    
    # Fit the model
    fit <- relmatLmer(value ~ (1|snp), df1, relmat = list(snp = K))
    
    # Heritability
    herit_boot(object = fit, exp = "snp / (snp + Residual)", boot.reps = 100)
  })

## Plot
g_fw_herit <- S2_MET_fw_herit %>%
  ggplot(aes(x = trait, y = heritability - bias, fill = term, group = term)) +
  geom_col(position = "dodge") +
  geom_col(aes(y = heritability), fill = "grey") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), position = position_dodge(0.9),
                width = 0.5) +
  theme_bw()



