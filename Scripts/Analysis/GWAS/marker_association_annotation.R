## Annotate GWAS results
## 
## This script will include code for doing the following:
## 1. Find overlap between associations for the mean and stability
## 2. Find overlap between associations detected in our analysis versus
## results from previous analyses
## 3. Look at the LD pattern surrounding significant associations
## 4. Generate examples of overlap and nearby gene annotations.
## 

# Load some packages
library(tidyverse)
library(readxl)
library(GenomicRanges)

# Directory containing the project repository
repo_dir <- getwd()

# Project and other directories
source(file.path(repo_dir, "source.R"))

# Load the genome annotation intersections
load(file.path(result_dir, "snp_genome_annotation.RData"))


# Read in the QTL metadata files that were downloaded from T3
qtl_meta <- map_df(list.files(data_dir, pattern = "meta", full.names = TRUE), read_csv)

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

# Read in the significant markers identified in our analysis
load(file.path(result_dir, "gwas_mlmm_marker_associations.RData"))

## Create GRanges for each marker data.frame
qtl_meta_grange <- qtl_meta_use %>%
  split(.$trait) %>%
  map(~makeGRangesFromDataFrame(., keep.extra.columns = TRUE, start.field = "position", 
                                end.field = "position"))





# Add a window of 5 Mbp to either side of the significant markers
window <- 5000000

gwas_mlmm_grange <- gwas_mlmm_marker %>%
  mutate(start = pos - window, end = pos + window) %>%
  split(.$trait) %>%
  map(~makeGRangesFromDataFrame(., keep.extra.columns = TRUE))

## For each trait, see how many of the markers detected in our analysis overlap
## with those on T3
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
      select(trait, coef, marker, chrom = seqnames, pos, start, end, MAF:WA, 
             qtl_marker = marker1, qtl_pos = start1, gene, feature) %>% 
      arrange(trait, coef, chrom, pos)

    
  })

## Write to a table
save_file <- file.path(fig_dir, "gwas_significant_associations.csv")
write_csv(x = gwas_markers_overlap, path = save_file)


# Read in the annotated excel file
gwas_markers_overlap <- read_excel(path = file.path(fig_dir, "Main_Figures/gwas_significant_associations.xlsx"),
                                   na = c("NA")) %>%
  mutate(qtl_pos = parse_number(qtl_pos))

# Number of associations per trait



## Count the number of overlaps per significant GWAS marker
gwas_markers_overlap_count <- gwas_markers_overlap %>% 
  group_by(trait, coef, marker) %>% 
  summarize(n_overlap = sum(!is.na(source)))

## Summarize overlaps per trait-coef combination
gwas_markers_overlap_count %>% 
  summarize(count_overlap = sum(n_overlap > 0),
            prop_overlap = mean(n_overlap > 0))

gwas_markers_overlap_count %>% 
  ungroup() %>% 
  filter(coef != "Genotype Mean") %>%   
  summarize(count_overlap = sum(n_overlap > 0),
            prop_overlap = mean(n_overlap > 0))


# What proportion of markers associated with stability overlapped with previous
# main-effect QTL?
# 

# Average distance of overlapping markers with the closest overlapping QTL
gwas_markers_overlap %>% 
  # Calculate the distance from the GWAS marker to the QTL marker
  mutate(marker_qtl_distance = abs(pos - qtl_pos)) %>%
  group_by(trait, coef, marker) %>% 
  summarize(marker_qtl_distance = min(marker_qtl_distance)) %>%
  summarize(mean_distance = mean(marker_qtl_distance, na.rm = TRUE))








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
  theme_manhattan(),
  theme(panel.border = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.line.y = element_line())
)


### Example for heading date

# Extract the results for mean per se and linear stability
gwas_hd_example <- gwas_pheno_mean_fw_adj %>% 
  filter(trait == "HeadingDate", coef %in% c("g", "b"), chrom == 5, between(pos, 575000000, 625000000))

# First plot the results for mean per se
g_hd_example_mean <- gwas_hd_example %>% 
  filter(coef == "g") %>% 
  ggplot(aes(x = pos / 1000000, y = -log10(q_value))) + 
  g_mod_man +
  ylab(expression(log[10](italic(q))~'for mean'~italic(per~se)))

# Second plot the results for linear stability
g_hd_example_fw <- gwas_hd_example %>% 
  filter(coef == "b") %>% ggplot(aes(x = pos / 1000000, y = log10(q_value))) + 
  g_mod_man +
  ylab(expression(log[10](italic(q))~'for linear stability'))

# Plot gene annotations
gene_annotations <- snp_info_grange_overlap$genes %>% 
  as.data.frame() %>% 
  select(chrom = 1, gene_start = ..start, gene_end = ..end, gene_id) %>% 
  mutate(chrom = parse_number(chrom)) %>% 
  filter(chrom == 5, between(gene_start, 575000000, 625000000), 
         between(gene_end, 575000000, 625000000)) %>% 
  distinct()

# Data.frame for Vrn1
vrn1_data <- data.frame(chrom = 5, gene_start = 599135017, gene_end = 599147377, gene_id = "Vrn-H1")

g_hd_annotations <- gene_annotations %>% 
  ggplot(aes(x = gene_start / 1000000, xend = gene_end / 1000000, y = 0, yend = 0)) + 
  geom_segment(lwd = 5) + 
  # Add VRN-H1
  geom_segment(data = vrn1_data, 
               aes(x = (gene_start - 100000) / 1000000, xend = (gene_end + 100000) / 1000000,
                   y = 0.1, yend = 0.1), inherit.aes = T, lwd = 5) + 
  geom_label(data = vrn1_data, aes(label = gene_id, y = 0.1), inherit.aes = T, vjust = 1.5, fontface = "italic") +
  xlab("Position on Chromosome 5 (Mbp)") +
  theme_bw() +
  theme(axis.line.x = element_line(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank())

# 
# ## Subset the pleiotropy results
# g_hd_example_plei <- gwas_pheno_mean_fw_plei_toplot %>% 
#   filter(trait == "HeadingDate", stab_coef == "Linear Stability", chrom == 5, 
#          between(pos, 575000000, 625000000)) %>%
#   ggplot(aes(x = pos / 1000000, y = -log10(qvalue), color = color)) + 
#   xlab("Position on Chromosome 5 (Mbp)") +
#   g_mod_man +
#   theme(axis.text.x = element_text(),
#         axis.ticks.x = element_line(),
#         axis.line.x = element_line(),
#         axis.title.x = element_text())


# Combine
hd_example_combine <- plot_grid(g_hd_example_mean, g_hd_annotations, g_hd_example_fw, 
                                ncol = 1, align = "hv", 
                                rel_heights = c(0.5, 0.2, 0.5))

## Save
save_file <- file.path(fig_dir, "hd_annotation_example.jpg")
ggsave(filename = save_file, plot = hd_example_combine, width = 6, height = 8, dpi = 1000)



## Plot the example for height
tr <- "PlantHeight"
chr <- 7
start <- 615e6
end <- 625e6

# Extract the results for mean per se and linear stability
gwas_ph_example <- gwas_pheno_mean_fw_adj %>% 
  filter(trait == tr, coef %in% c("g", "b"), chrom == chr, between(pos, start, end))

# First plot the results for mean per se
g_ph_example_mean <- gwas_ph_example %>% 
  filter(coef == "g") %>% 
  ggplot(aes(x = pos / 1000000, y = -log10(q_value))) + 
  g_mod_man +
  ylab(expression(log[10](italic(q))~'for mean'~italic(per~se)))

# Second plot the results for linear stability
g_ph_example_fw <- gwas_ph_example %>% 
  filter(coef == "b") %>% ggplot(aes(x = pos / 1000000, y = log10(q_value))) + 
  g_mod_man +
  ylab(expression(log[10](italic(q))~'for linear stability'))

gene_annotations <- snp_info_grange_overlap$genes %>% 
  as.data.frame() %>% 
  select(chrom = 1, gene_start = ..start, gene_end = ..end, gene_id) %>% 
  mutate(chrom = parse_number(chrom)) %>% 
  filter(chrom == chr, between(gene_start, start, end), between(gene_end, start, end)) %>% 
  distinct()

# Add annotation
g_ph_annotations <- gene_annotations %>% 
  ggplot(aes(x = gene_start / 1000000, xend = gene_end / 1000000, y = 0, yend = 0)) + 
  geom_segment(lwd = 5) + 
  # # Add VRN-H1
  # geom_segment(data = vrn1_data, 
  #              aes(x = (gene_start - 100000) / 1000000, xend = (gene_end + 100000) / 1000000,
  #                  y = 0.1, yend = 0.1), inherit.aes = T, lwd = 5) + 
  # geom_label(data = vrn1_data, aes(label = gene_id, y = 0.1), inherit.aes = T, vjust = 1.5, fontface = "italic") +
  xlab("Position on Chromosome 5 (Mbp)") +
  theme_bw() +
  theme(axis.line.x = element_line(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),
        panel.grid = element_blank())

# Combine
ph_example_combine <- plot_grid(g_ph_example_mean, g_ph_annotations, g_ph_example_fw, 
                                ncol = 1, align = "hv", 
                                rel_heights = c(0.5, 0.2, 0.5))










#### Intervals around significant associations

# Group significant SNPs for each trait that are within 5 Mbp and in high LD ($r^2 \ge 0.80$). 


## Write a function that uses distance and LD to group SNPs
snp_clust <- function(grange, LD, wind.size, LD.threshold) {
  
  # Grab the SNP positions
  pos <- start(grange)
  # Grab the marker names
  mar_names <- names(grange)
  
  # Start some vectors
  dist_list <- LD_list <- vector("numeric", length(grange))
  
  # Start a distance group and LD group count
  dist_group <- LD_group <- 1
  
  # Add the first SNP to group 1
  dist_list[1] <- LD_list[1] <- dist_group
  
  # Start the snp counter
  i = 2
  # Iterate over the markers in order
  while (i <= length(grange)) {
    
    # Is marker i within the window size of the previous marker?
    window <- c(pos[i] - wind.size, pos[i] + wind.size)
    is_within <- between(pos[i - 1], left = window[1], right = window[2])
    
    # If so, add this marker to the dist group of the previous
    if (is_within) {
      dist_list[i] <- dist_group
    } else {
      dist_group <- dist_group + 1
      dist_list[i] <- dist_group
    }
    
    i_dist_group <- setdiff(which(dist_list == dist_group), i)
    # Is marker i within the LD cutoff of any of the markers in the same distance group?
    is_in_LD <- any(LD[i, i_dist_group] >= LD.threshold)
    
    # If so, add this marker to the LD group of the previous
    if (is_in_LD) {
      LD_list[i] <- LD_group
    } else {
      LD_group <- LD_group + 1
      LD_list[i] <- LD_group
    }
    
    # enumerate i
    i <- i + 1
    
  }
  
  # Group
  data.frame(marker = mar_names, group = as.numeric(as.factor(dist_list + LD_list)), 
             stringsAsFactors = FALSE)
  
}

# LD threshold for grouping markers
LD_threshold <- 0.80
# Window size for genomic region
# Note that half of this value will be used for overlaps
wind_size <- 10000000

# Create GRanges objects
sig_markers_GRange <- gwas_mlmm_Gmodel %>%
  ungroup() %>%
  mutate(snp_pos = pos) %>%
  group_by(trait, coef, chrom) %>%
  do(grange = {
    df <- . 
    makeGRangesFromDataFrame(
      df = df, keep.extra.columns = TRUE, ignore.strand = TRUE, seqnames.field = "chrom", 
      start.field = "pos", end.field = "pos") %>% 
      `ranges<-`(., value = IRanges(start = start(.) - (wind_size / 2), end = end(.) + (wind_size / 2))) })

# Calculate LD between all markers for each group,
# then use the LD to cluster markers
sig_markers_groups <- sig_markers_GRange %>% 
  group_by(trait, coef, chrom) %>%
  do(grange = {
    grange <- .$grange[[1]]
    
    # If only one SNP, return a regular data.frame
    if (length(grange) == 1) {
      grange$group <- "group1"
      
    } else {
      
      # Calculate LD
      LD <- cor(M[,grange$marker])^2
      
      # Add names to the grange
      names(grange) <- grange$marker
      
      # Cluster
      groups <- snp_clust(grange = grange, LD = LD, wind.size = wind_size, LD.threshold = LD_threshold)
      grange$group <- str_c("group", groups$group)
    }
    
    grange
  })

# Combine the data
sig_markers_grange_df <- sig_markers_groups$grange %>% 
  map_df(as.data.frame) %>% 
  as_data_frame() %>%
  dplyr::rename(chrom = seqnames)

# Save this
save_file <- file.path(result_dir, "S2MET_gwas_sig_loci.RData")
save("sig_markers_grange_df", file = save_file)






# Summarize the number of significant markers and intervals

## Summarize the number of significant markers and the number of sgnificant intervals
sig_summary <- sig_markers_grange_df %>% 
  group_by(trait, coef, chrom) %>% 
  summarize(SNP = n_distinct(marker), Interval = n_distinct(group)) %>% 
  summarize_at(vars(SNP, Interval), sum)

# Translator for coefficient
coef_repl <- function(x) {
  case_when(x == "log_delta" ~ "Non-Linear Stability",
            x == "g" ~ "Genotype Mean", 
            x == "b" ~ "Linear Stability")
}

## Plot the intervals
# First change the coefficient to a factor and create y axis 
g_sig_interval <- sig_markers_grange_df %>% 
  mutate(coef = parse_factor(coef, levels = c("g", "b", "log_delta")),
         y_coef = as.numeric(coef) / 2) %>% 
  ggplot() + 
  geom_segment(aes(x = y_coef, xend = y_coef, y = start / 1000000, 
                   yend = end / 1000000, col = coef), size = 3) + 
  geom_segment(data = barley_lengths, aes(x = 0, y = 0, xend = 0, yend = length / 1000000), 
               inherit.aes = FALSE, size = 2, lineend = "round", col = "grey70") + 
  ylab("Position (Mbp)") +
  xlim(c(-1, 2)) +
  facet_grid(trait ~ chrom, switch = "x") + 
  scale_y_reverse() +
  theme_bw() +
  theme(panel.grid = element_blank(), 
        panel.border = element_blank(), 
        panel.spacing.y = unit(0.5, units = "cm"),
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(),
        strip.background = element_blank())

save_file <- file.path(fig_dir, "significant_loci_intervals.jpg")
ggsave(filename = save_file, plot = g_sig_interval, height = 7, width = 9)



# Create genomic ranges for each of the trait-coef combinations

grange_sep <- sig_markers_grange_df %>% 
  group_by(trait, coef) %>% 
  do(Grange = GRanges(.)) %>%
  ungroup()

# Group by trait and find potential overlap between genotype mean intervals and
# plasticity intervals
grange_compare <- grange_sep %>% 
  group_by(trait) %>% 
  do({
    df <- .
    
    # Data.frame of combinations
    combn_df <- combn(x = df$coef, m = 2) %>%
      t() %>%
      as.data.frame() %>%
      `colnames<-`(., c("coef1", "coef2")) %>%
      map(~as.character(.))
    
    # Iterate over the combinations
    compare_out <- map2(.x = combn_df$coef1, .y = combn_df$coef2, 
                        function(c1, c2) {
                          grange1 <- subset(df, coef == c1)$Grange[[1]]
                          grange2 <- subset(df, coef == c2)$Grange[[1]]
                          
                          # Calculate the distance to neareast and add the grange information
                          dtn <- suppressWarnings(distanceToNearest(x = grange1, subject = grange2)) %>%
                            as_data_frame() %>% 
                            mutate(queryGrange = map(queryHits, ~grange1[.]), 
                                   subjectGrange = map(subjectHits, ~grange2[.]))
                        })
    
    # Return the comparisons
    combn_df %>% 
      as_data_frame() %>%
      mutate(comparison = compare_out) })


# Look at the comparisons where overlap is present and determine the LD between
# the SNP in those intervals
grange_compare_LD <- grange_compare %>%
  unnest() %>% 
  filter(distance == 0) %>%
  do({
    
    df <- .
    
    # Extract the markers
    overlap_df_grange <- df %>% 
      select(trait:coef2, contains("Grange")) %>%
      gather(temp, coef, coef1, coef2) %>% 
      gather(temp1, Grange, contains("Grange")) %>% 
      filter(temp == "coef1" & temp1 == "queryGrange" | temp == "coef2" & temp1 == "subjectGrange") %>% 
      select(-contains("temp")) %>% 
      mutate(Grange = map(Grange, as.data.frame)) %>% 
      unnest()
    
    # Return the LD
    overlap_df_grange %>% 
      select(trait, coef, chrom = seqnames, start, end, marker:group) %>% 
      mutate(LD = cor(M[,overlap_df_grange$marker])[1,2]^2)
    
  })

grange_compare_LD %>% 
  select(trait, coef, chrom, alpha, q_value, snp_pos, LD)




## For each of the significant intervals, plot the LD of those markers with nearby markers within 10 Mbp

# In each group, select the marker with the highest p-value
sig_markers_grange_df_select <- sig_markers_grange_df %>% 
  group_by(trait, coef, chrom, group) %>% 
  filter(p_value == min(p_value))

# Calculate LD across all markers on each chromosome
snp_LD <- snp_info %>% 
  split(.$chrom) %>% 
  map("marker") %>% 
  map(~cor(M[,.])^2)

# Iterate over each marker and calculate the LD between that marker
# and all nearby markers within the interval
sig_markers_grange_df_select_LD <- sig_markers_grange_df_select %>%
  group_by(trait, coef, marker) %>%
  do({
    df <- .
    
    nearby_markers <- snp_info %>% 
      filter(chrom == df$chrom, between(pos, left = df$start, right = df$end),
             marker != df$marker)
    
    # Extract LD
    df_LD <- snp_LD[[df$chrom]][df$marker, nearby_markers$marker]
    
    # Add the LD to the nearby markers and return
    to_return <- nearby_markers %>% 
      mutate(LD = df_LD) 
    
    df %>% 
      select(chrom, snp_pos, alpha:q_value) %>% 
      mutate(snp_LD = list(to_return)) })

# Plot the significant marker and the LD of the surrounding markers
g_sig_marker_LD <- sig_markers_grange_df_select_LD %>% 
  unnest() %>% 
  ggplot() + 
  geom_point(aes(x = pos / 1000000, y = LD)) + 
  geom_vline(aes(xintercept = snp_pos / 1000000)) + 
  ylab(expression(r^2)) + 
  xlab("Position (Mbp)") + 
  ylim(c(0,1)) +
  facet_wrap_paginate(~ trait + coef + marker, nrow = 3, ncol = 3, scales = "free_x")

# Iterate over the pages and save







## Window overlap
## 
## First find the overlap between mean and stability and quantify
## Then find the overlap between stability and marker effect stability

# Convert the mean and stability intervals into GRanges
sig_intervals <- sig_loci_df_intervals %>% 
  droplevels() %>% 
  mutate(seqnames = str_c("chr", .$chrom, "H")) %>%
  distinct(trait, assoc_test, interval_start, interval_end, seqnames) %>%
  group_by(assoc_test, trait) %>% 
  do(range = GRanges(seqnames = .$seqnames, ranges = IRanges(start = .$interval_start, end = .$interval_end))) %>%
  select(trait, names(.)) %>%
  spread(assoc_test, range)



# Write a function that takes two lists of Ranges and finds the overlaps
# and number of overlaps
find_overlaps <- function(q, s) {
  # If null, return NA
  if (is.null(q) | is.null(s)) {
    data_frame(hits = list(NA), query_hits = list(NA), subject_hits = list(NA),
               query_length = as.integer(NA), subject_length = as.integer(NA),
               query_hits_n = as.integer(NA), subject_hits_n = as.integer(NA))
    
    
  } else {
    # Find overlaps
    hits <- suppressWarnings(findOverlaps(query = q, subject = s))
    # Data frame of query and subject hits, lengths, number
    data_frame(hits = list(hits)) %>% 
      mutate(query_hits = map(hits, ~q[queryHits(.),]), 
             subject_hits = map(hits, ~s[subjectHits(.),]),
             query_length = length(q), subject_length = length(s)) %>% 
      mutate_at(vars(contains("_hits")), funs(n = map_dbl(., length)))
  }
}



# For each trait, find the overlaps between mean and stability
interval_overlaps <- sig_intervals %>% 
  mutate(mean_linear_stability_overlap = list(main_effect_mean, linear_stability) %>% pmap(find_overlaps),
         mean_nonlinear_stability_overlap = list(main_effect_mean, `non-linear_stability`) %>% pmap(find_overlaps),
         linear_stability_mar_linear_stability_overlap = list(linear_stability, linear_marker_stability) %>% pmap(find_overlaps),
         linear_stability_mar_linear_plasticity_overlap = list(linear_stability, linear_marker_plasticity) %>% pmap(find_overlaps),
         nonlinear_stability_mar_nonlinear_stability_overlap = list(`non-linear_stability`, `non-linear_marker_plasticity`) %>% pmap(find_overlaps) )

# Extract the data.frame results
interval_overlaps_df <- interval_overlaps %>% 
  select(trait, contains("overlap")) %>% 
  gather(overlap, data, -trait) %>% 
  unnest() %>% 
  # Calculate proportions
  mutate(prop_query_hits = query_hits_n / query_length,
         prop_subject_hits = subject_hits_n / subject_length)

interval_overlaps_df %>% select(trait, overlap, query_hits_n:prop_subject_hits)






