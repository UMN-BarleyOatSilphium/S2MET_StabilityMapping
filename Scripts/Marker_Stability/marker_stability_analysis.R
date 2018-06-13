## S2MET Mapping
## Calculate marker effect reaction norms
## 
## 
## Author: Jeff Neyhart
## Last updated: June 12, 2018
## 


# Load the source
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

# Load the stability results
load(file.path(result_dir, "pheno_mean_fw_results.RData"))
# Load the marker effect by environment results
load(file.path(result_dir, "marker_by_env_effects.RData"))


# Get the marker information
snp_info <- S2TP_imputed_multi_genos_hmp %>%
  select(marker = rs, chrom, pos, cM_pos) %>%
  # Correct the cM position for BOPA snps
  mutate(cM_pos = if_else(str_detect(marker, "^S"), cM_pos, cM_pos / 1000))


M <- S2TP_imputed_multi_genos_mat




## Manage the mxe data and convert to a tidy data.frame
mxe_df <- marker_by_env_effects %>% 
  list(., names(.)) %>% 
  pmap_df(~mutate(.x, trait = .y)) %>%
  select(trait, environment, names(.))

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
marker_mean_fw <- mxe_df1 %>%
  full_join(., marker_effect_stab, by = c("trait", "marker")) %>%
  left_join(., snp_info, by = "marker") %>%
  select(marker, chrom:cM_pos, trait, environment, effect, h:delta, df)

# Save the data
save_file <- file.path(result_dir, "S2MET_marker_mean_fw_results.RData")
save("marker_mean_fw", "marker_effect_stab", file = save_file)





###### Below here, data from above can be loaded

load(file.path(result_dir, "S2MET_marker_mean_fw_results.RData"))


## Plot the distribution of stability estimates

## Transform
marker_mean_fw_trans <- marker_mean_fw %>% 
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
g_marker_fw_dens_b <- marker_mean_fw_trans %>%
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
# g_marker_fw_dist <- g_marker_fw_dens_b + g_marker_fw_dens_delta

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
xtabs(~ trait + significance, S2_MET_marker_eff_pheno_fw_sig)




## Plot the stability estimates across the genome with the empirical threshold
# Color theme
colors <- c("black", umn_palette(n = 4)[-c(1:2)]) %>%
  set_names(c("Average", "Plastic", "Stable"))


g_mar_stab <- S2_MET_marker_eff_pheno_fw_sig %>%   
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


## Poster version
g_mar_stab_poster <- g_mar_stab +
  geom_point(size = 2) + 
  labs(caption = "95% empirical threshold given by dotted line.") +
  theme_poster() + 
  theme(legend.position = "bottom",
        panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.text = element_text(size = 14),
        panel.spacing.x = unit(x = 0, units = "cm"),
        panel.spacing.y = unit(x = 1, units = "lines"),
        panel.border = element_rect(color = "grey75"),
        plot.caption = element_text(size = 10))

ggsave(filename = "marker_stability_genome_poster.jpg", plot = g_mar_stab_poster, path = fig_dir,
       height = 10, width = 16, dpi = 1000)





## Gene Annotation

# Look at gene annotation with relation to the different groups of SNPs. 
# Determine the proportion of SNPs in genes, gene-proximal, or intergenic


## Load the barley annotation
ann_dir <- "C:/Users/Jeff/GoogleDrive/BarleyLab/Projects/Genomics/Annotation"
load(file.path(ann_dir, "barley_genomic_ranges.RData"))


## First plot gene annotation density across the genome
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
g_plast_gene_dens <- plot_grid(g_mar_stab_poster + theme(axis.title.x = element_blank()), 
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









#### Overlap with Stability and Mean Loci
## Also look at previously described QTL

# Load the GWAS results for genotypic mean and phenotypic stability and see if the 
# markers with the highest plasticity/stability measure overlap
 
load(file.path(result_dir, "gwas_adjusted_significant_results.RData"))

# Read in the QTL metadata files that were downloaded from T3
qtl_meta <- map_df(list.files(data_dir, pattern = "meta", full.names = TRUE), read_csv)
# Read in association data from Pauli2014 and Wang2012
bcap_association <- read_csv(file.path(data_dir, "BCAP_association_qtl.csv")) %>%
  mutate(position = parse_number(position)) %>%
  select(trait, marker, chrom = chromosome, pos = position, gene:feature, reference) %>%
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
  select(trait, marker, chrom, pos = position, gene, feature) %>%
  arrange(trait, chrom, pos)

# Combine data
qtl_meta_use1 <- bind_rows(qtl_meta_use, bcap_association)



### First highlight some examples for known genes / QTL for heading date and plant height

## Heading Date
tr <- "HeadingDate"
chr <- 5
start <- 580e6
end <- 620e6


## Add information on VRN
vrn1_data <- data.frame(chrom = 5, gene_start = 599135017, gene_end = 599147377, gene_id = "Vrn-H1")

# Get the QTL information
hd_qtl <- qtl_meta_use1 %>% 
  filter(trait == tr, chrom == chr)



## Subset the chromosome of interest
g_hd_marstab <- S2_MET_marker_eff_pheno_fw_sig %>%   
  filter(coef == "b", trait == tr, chrom == chr, between(pos, start, end)) %>%
  ggplot(aes(x = pos / 1000000, y = estimate, col = significance)) + 
  geom_point() + 
  geom_hline(aes(yintercept = lower_perc), lty = 2) +
  geom_hline(aes(yintercept = upper_perc), lty = 2) +
  scale_color_manual(values = colors, name = NULL) +
  ylab("Marker Effect\nPlasticity Estimate") +
  xlab("Position (Mbp)") +
  labs(title = "Heading Date") + 
  theme_poster() +
  theme(legend.position = c(0.15, 0.1),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_blank())


## Create a separate gene annotation plot
g_hd_ann <- vrn1_data %>%
  ggplot(aes(y = 1, yend = 1)) +
  xlim(c(start / 1e6, end / 1e6)) + 
  geom_segment(aes(x = (gene_start - 100000) / 1000000, xend = (gene_end + 100000) / 1000000), lwd = 5) +
  geom_text(aes(x = gene_start / 1000000, label = gene_id), vjust = 2, fontface = "italic") +
  # Add significant QTL
  geom_segment(data = hd_qtl, aes(x = (pos - 100000) / 1000000, xend = (pos + 100000) / 1000000), lwd = 5) +
  xlab("Position on Chromosome 5 (Mbp)") +
  ylab("Known\nGenes/QTL") +
  ylim(c(0, 1)) +
  theme_poster() +
  theme(axis.line.x = element_line(),
        axis.text.y = element_blank(),
        # axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank())

## Add the GWAS analysis of stabililty
g_hd_stab_gwas <- gwas_pheno_mean_fw_adj %>% 
  filter(trait == tr, coef == "b", chrom == chr, between(pos, start, end)) %>%
  ggplot(aes(x = pos / 1000000, y = qvalue_neg_log10)) + 
  geom_point() +
  ylab(expression(atop(-log[10](italic(q)),'Phenotypic Plasticty'))) +
  theme_poster() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_blank())
  
  
## Combine
g_hd_example <- plot_grid(
  g_hd_marstab, 
  g_hd_stab_gwas,
  g_hd_ann,
  ncol = 1, rel_heights = c(0.8, 0.4, 0.25), align = "hv", axis = "lr")

## Save
ggsave(filename = "hd_marstab_annotation_poster.jpg", plot = g_hd_example, path = fig_dir,
       height = 8, width = 6, dpi = 1000)


## Plant Height
tr <- "PlantHeight"
chr <- 6
start <- 0
end <- 20e6


ph_qtl <- qtl_meta_use1 %>% 
  filter(trait == tr, chrom == chr)

## Plot of marker effect plasticity
g_ph_marstab <- S2_MET_marker_eff_pheno_fw_sig %>%   
  filter(coef == "b", trait == tr, chrom == chr, between(pos, start, end)) %>%
  ggplot(aes(x = pos / 1000000, y = estimate, col = significance)) + 
  geom_point() + 
  geom_hline(aes(yintercept = lower_perc), lty = 2) +
  geom_hline(aes(yintercept = upper_perc), lty = 2) +
  scale_color_manual(values = colors, name = NULL) +
  ylab("Marker Plasticity Estimate") +
  xlab("Position (Mbp)") +
  theme_poster() +
  theme(legend.position = c(0.15, 0.05),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_blank())


## Create a separate gene annotation plot
g_ph_ann <- ph_qtl %>%
  ggplot(aes(y = 0, yend = 0)) +
  xlim(c(start / 1e6, end / 1e6)) + 
  geom_segment(aes(x = (pos - 100000) / 1000000, xend = (pos + 100000) / 1000000), lwd = 5) +
  xlab("Position on Chromosome 6 (Mbp)") +
  theme_poster() +
  theme(axis.line.x = element_line(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank())

## Add the GWAS analysis of stabililty
g_ph_stab_gwas <- gwas_pheno_mean_fw_adj %>% 
  filter(trait == tr, coef == "b", chrom == chr, between(pos, start, end)) %>%
  ggplot(aes(x = pos / 1000000, y = qvalue_neg_log10)) + 
  geom_point() +
  ylab(expression(-log[10](italic(q))~'Phenotypic Plasticty')) +
  theme_poster() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.border = element_blank())

## Combine
g_ph_example <- plot_grid(
  g_ph_marstab, 
  g_ph_stab_gwas,
  g_ph_ann,
  ncol = 1, rel_heights = c(0.8, 0.4, 0.2), align = "hv", axis = "lr")

## Save
ggsave(filename = "ph_marstab_annotation_poster.jpg", plot = g_ph_example, path = fig_dir,
       height = 8, width = 5, dpi = 1000)


## Combine HD and PH
g_hd_ph_example <- plot_grid(g_hd_example, g_ph_example, ncol = 2)
## Save
ggsave(filename = "hd_ph_marstab_annotation_poster.jpg", plot = g_hd_ph_example, path = fig_dir,
       height = 8, width = 10, dpi = 1000)





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
K <- A.mat(X = S2TP_imputed_multi_genos_mat, min.MAF = 0, max.missing = 1)

# Estimate the proportion of variation due to all SNPs
S2_MET_mean_fw_varcomp <- S2_MET_pheno_fw_uniq_trans %>% 
  group_by(trait, coef) %>%
  do({
    # Extract data
    df <- .
    
    # Copy the line_name column
    df1 <- df %>% mutate(snp = line_name)
    
    # Fit the model
    fit <- relmatLmer(value ~ (1|snp), df1, relmat = list(snp = K))
    
    # Extract the proportion of variance attributed to snps and to residuals
    VarProp(fit) %>% 
      select(-contains("var"))
  })

# Remove extra data
mean_fw_varcomp <- S2_MET_mean_fw_varcomp %>%
  filter(grp == "snp") %>%
  select(trait, coef, prop) %>%
  ungroup()


# Number of model fittings
n_iter <- 100

# Now for each trait, find the proportion of variation in the mean and stability
# due to different marker groups
S2_MET_mean_fw_martype_varcomp <- S2_MET_pheno_fw_uniq_trans %>%
  group_by(trait, coef) %>%
  do({
    
    df <- .
    
    # What is the trait we are dealing with
    tr <- unique(df$trait)
    
    # Extract the markers for the particular trait and separate by average/stable/plastic
    marker_types <- S2_MET_marker_eff_pheno_fw_sig %>%
      filter(trait == tr) %>% 
      split(.$significance) %>% 
      map("marker")
    
    # Use the plastic and stable markers to make relationship matrices
    K_plas <- marker_types$Plastic %>% M[,.,drop = F] %>% A.mat(X = ., min.MAF = 0, max.missing = 1)
    K_stab <- marker_types$Stable %>% M[,.,drop = F] %>% A.mat(X = ., min.MAF = 0, max.missing = 1)
    
    # Sample size
    sample_size <- length(marker_types$Plastic)
    
    # Create random samples of the average markers
    average_marker_samples <- rerun(.n = n_iter, sample(marker_types$Average, size = sample_size))
    
    ## Fit using sommer
    y <- df$value
    X <- model.matrix(~ 1, df)
    Z_g <- model.matrix(~ -1 + line_name, df) %>%
      `colnames<-`(., colnames(K_stab))
    
    # Iterate over the samples
    var_comp_out <- average_marker_samples %>%
      map(function(marker_sample) {
        
        # Create a relationship matrix
        K_avg <- M[,marker_sample,drop = F] %>% A.mat(X = ., min.MAF = 0, max.missing = 1)
        
        # Create a Z list
        Z <- list(
          stab = list(Z = Z_g, K = K_stab),
          plas = list(Z = Z_g, K = K_plas),
          avg = list(Z = Z_g, K = K_avg)
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


## The results from mapping suggest that stable and sensitive markers reside in 
## similar locations (or some overlap exists). This suggest that these markers may
## tag the same region. If this is true, we would expect a combination of stable
## and plastic markers to not explain much more of the variation as either
## one. Below we test that

S2_MET_mean_fw_comb_martype_varcomp <- S2_MET_pheno_fw_uniq_trans %>%
  group_by(trait, coef) %>%
  do({
    
    df <- .
    
    # What is the trait we are dealing with
    tr <- unique(df$trait)
    
    # Extract the markers for the particular trait and separate by average/stable/plastic
    marker_types <- S2_MET_marker_eff_pheno_fw_sig %>%
      filter(trait == tr) %>% 
      split(.$significance) %>% 
      map("marker")
    
    # Sample size
    sample_size <- length(marker_types$Plastic)
    
    # Create random samples of the average markers
    average_marker_samples <- rerun(.n = n_iter, sample(marker_types$Average, size = sample_size))
    # Create random samples of a 50-50 mix of stable and plastic markers
    mix_marker_samples <- rerun(.n = n_iter, c(sample(marker_types$Plastic, size = sample_size / 2), 
                                               sample(marker_types$Stable, size = sample_size / 2)))
    
    ## Fit using sommer
    mf <- model.frame(value ~ line_name, df)
    y <- model.response(mf)
    X <- model.matrix(~ 1, mf)
    Z_g <- model.matrix(~ -1 + line_name, mf) %>%
      `colnames<-`(., mf$line_name)
    
    # Iterate over the samples
    var_comp_out <- list(average_marker_samples, mix_marker_samples) %>%
      pmap(function(avg_sample, mix_sample) {
        
        # Create a relationship matrix
        K_avg <- M[,avg_sample,drop = F] %>% A.mat(X = ., min.MAF = 0, max.missing = 1)
        # Create mixed relationship matrix
        K_mix <- M[,mix_sample,drop = F] %>% A.mat(X = ., min.MAF = 0, max.missing = 1)
        
        # Create a Z list
        Z <- list(
          mix = list(Z = Z_g, K = K_mix),
          avg = list(Z = Z_g, K = K_avg)
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



## Character vector for replacing marker types
mar_type_replace <- c("avg" = "Average", "plas" = "Plastic", "stab" = "Stable")


## Summarize
# Separate variance components for stable, plastic, and average markers
S2_MET_mean_fw_martype_varcomp_summ <- S2_MET_mean_fw_martype_varcomp %>% 
  mutate(total = stab + plas + avg + units,
         stab_plas = stab + plas) %>% # Add stab + plas together
  mutate_at(vars(stab:units, stab_plas), funs(. / total)) %>%
  select(-units, -total) %>% 
  gather(marker_type, prop_var, stab:avg, stab_plas) %>%
  ungroup()

## Plot
S2_MET_mean_fw_martype_varcomp_summ %>% 
  filter(marker_type %in% names(mar_type_replace)) %>%
  mutate(coef = str_replace_all(coef, rev(coef_replace)),
         coef = factor(coef, levels = coef_replace),
         marker_type = str_replace_all(marker_type, mar_type_replace)) %>%
  ggplot(aes(x = coef, y = prop_var, col = marker_type)) + 
  geom_boxplot(position = "dodge")




# Separate variance components for 50/50 stable/plastic markers and average markers
S2_MET_mean_fw_comb_martype_varcomp_summ <- S2_MET_mean_fw_comb_martype_varcomp %>% 
  mutate(total = mix + avg + units) %>% 
  mutate_at(vars(mix:units), funs(. / total)) %>%
  select(-units, -total) %>% 
  gather(marker_type, prop_var, mix:avg) %>%
  ungroup()


# Combine
mean_stability_martype_varcomp <- bind_rows(
  S2_MET_mean_fw_martype_varcomp_summ,
  filter(S2_MET_mean_fw_comb_martype_varcomp_summ, marker_type == "mix")) %>%
  # Convert factors
  mutate(marker_type = parse_factor(marker_type, levels = c("avg", "plas", "stab",  "stab_plas", "mix")))


# Group and take the mean and sd 
# Output in a table format            
mean_stability_martype_varcomp_table <- mean_stability_martype_varcomp %>% 
  group_by(trait, coef, marker_type) %>% 
  summarize_at(vars(prop_var), funs(mean, sd)) %>%
  mutate_at(vars(mean, sd), funs(round(., 3))) %>% 
  mutate(temp = str_c(mean, " (", sd, ")")) %>% 
  select(-mean, -sd) %>% 
  spread(marker_type, temp)

# Plot as a boxplot
mean_stability_martype_varcomp %>% 
  mutate(coef = str_replace_all(coef, c("b" = "Linear Stability", "log_delta" = "Non-Linear Stability", 
                                    "g" = "Genotype Mean")),
         coef = parse_factor(coef, levels = c("Linear Stability", "Non-Linear Stability", "Genotype Mean"))) %>%
  ggplot(aes(x = coef, y = prop_var, color = marker_type)) + 
  geom_boxplot(width = 1) + 
  facet_grid(trait ~ .) +
  theme_bw()




## Save
save_file <- file.path(result_dir, "S2MET_marker_stability_varcomp.RData")
save("mean_fw_varcomp", "mean_stability_martype_varcomp", file = save_file)




### Appendix
### 
### Plots and other code

## Plot examples of the stability of genotypes/markers
pheno_stable_ex <- S2_MET_pheno_mean_fw %>% 
  group_by(trait) %>% 
  filter(b == min(b) | b == max(b) | abs(b - 1) == min(abs(b - 1))) %>%
  mutate(class = case_when(
    b == max(b) ~ "Plastic",
    b == min(b) ~ "Stable",
    TRUE ~ "Average")) %>%
  ungroup()

g_pheno_stable_ex <- pheno_stable_ex %>% 
  # mutate(value = value - h) %>%
  ggplot(aes(x = h, y = value, color = class)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) +  
  scale_color_manual(values = colors, name = "Response Type") +
  facet_wrap(~ trait, ncol = 1, scales = "free") +
  ylab(expression("Phenotypic Value "~(italic(y[ij])))) + 
  xlab(expression("Environmental Effect "~(italic(t[j])))) +
  labs(title = "Phenotypic Plasticity") + 
  theme_poster() +
  theme(legend.position = "bottom", 
        legend.text = element_text(size = 14),
        title = element_text(size = 16))


mark_stable_ex <- S2_MET_marker_mean_fw %>%
  group_by(trait) %>% 
  filter(b == min(b) | b == max(b) | abs(b - 0) == min(abs(b - 0))) %>%
  mutate(class = case_when(
    b == max(b) ~ "Plastic",
    b == min(b) ~ "Stable",
    TRUE ~ "Average")) %>%
  ungroup()

g_mark_stable_ex <- mark_stable_ex %>% 
  # mutate(effect = effect + h) %>%
  ggplot(aes(x = h, y = effect, color = class)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) +  
  scale_color_manual(values = colors, name = NULL) +
  facet_wrap(~ trait, ncol = 1, scales = "free") +
  ylab(expression("Marker Effect "~(italic(alpha[jp])))) + 
  xlab(expression("Environmental Effect "~(italic(t[j])))) +
  labs(title = "Marker Effect Plasticity",
       caption = expression("Environmental effect"~(italic(t[j]))~"not added to marker effect to highlight reponse.")) + 
  theme_poster() +
  theme(legend.position = "none",
        plot.title = element_text(size = 16),
        plot.caption = element_text(size = 10))

## Combine
g_stable_ex <- plot_grid(g_pheno_stable_ex + theme(legend.position = "none"), 
                         g_mark_stable_ex, ncol = 2, align = "hv")
g_stable_ex1 <- plot_grid(g_stable_ex, get_legend(g_pheno_stable_ex), ncol = 1, rel_heights = c(0.95, 0.05))

## Save
ggsave(filename = "stability_example_poster.jpg", plot = g_stable_ex1, path = fig_dir,
       height = 10, width = 8, dpi = 1000)



















