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

