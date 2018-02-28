## LD pattern in the S2TP
## 
## This script will do two things:
## 1. Calculate and plot LD across the genome
## 2. Plot a heatmap of LD

# Set directories and such
library(tidyverse)
library(LDheatmap)
library(readxl)
library(broom)

# Project and other directories
source("C:/Users/Jeff/Google Drive/Barley Lab/Projects/S2MET_Mapping/source.R")

# Create a marker matrix M
M <- S2TP_imputed_multi_genos_mat

# Data frame of SNP information (pos, chrom, etc)
snp_info <- S2TP_imputed_multi_genos_hmp %>% 
  select(marker = rs, chrom, pos, cM_pos)

# Divide the SNPs by chromosome
snp_info1 <- snp_info %>% 
  split(.$chrom)

# Iterate over the chromosomes of SNPs and subset the marker matrix,
# then calculate pairwise r^2 and tidy it up
pairwise_LD <- snp_info1 %>% 
  map(~M[,.$marker] %>%
        {cor(.)^2} %>% 
        as.dist() %>% 
        tidy() %>% 
        rename(marker1 = item1, marker2 = item2, LD = distance) %>%
        as_data_frame() )

# Add the SNP information back in
pairwise_LD_snp_info <- list(pairwise_LD, snp_info1) %>%
  pmap(~left_join(x = .x, y = .y, by = c("marker1" = "marker")) %>%
         left_join(x = ., y = .y, by = c("marker2" = "marker")) %>%
         select(marker1:LD, chrom = chrom.x, pos1 = pos.x, pos2 = pos.y,
                cM_pos1 = cM_pos.x, cM_pos2 = cM_pos.y) ) %>%
  # Calculate the distance between the markers
  map(~mutate(., distance_bp = abs(pos1 - pos2), 
              distance_cM = abs(cM_pos1 - cM_pos2)))

# Using increments of 1 kbp, calculate the average LD between all pairs of 
# markers that are within that interval
intervals <- seq(from = 1000, to = 50000000, by = 1000)



# Calculate the average LD per interval per chromosome
pairwise_bpDist_LD <- pairwise_LD_snp_info %>%
  map(as.data.frame) %>%
  map_df(function(chrom) {
    # Empty vector to store results
    LD_vec <- numeric(length = length(intervals))
    # For loop
    for (i in seq_along(intervals)) {
      LD_vec[i] <- mean(subset(chrom, distance_bp <= intervals[i], LD, drop = TRUE)) }
    
    # Return data.frame
    data.frame(chrom = unique(chrom$chrom), interval = intervals, meanLD = LD_vec) })

# Filter all of the data for marker pairs that are at most max(intervals) apart
pairwise_LD_snp_info_plot <- pairwise_LD_snp_info %>%
  map_df(~filter(., distance_bp <= max(intervals)))


# Plot
pairwise_bpDist_LD %>%
  mutate(chrom = as.factor(chrom)) %>%
  ggplot() +
  geom_point(data = pairwise_LD_snp_info_plot, aes(x = distance_bp, y = LD), alpha = 0.5) +
  geom_point(aes(x = interval, y = meanLD, col = chrom))


## Caclculate pairwise LD as a function of the number of SNPs between pairs of SNPs
# First calculate the number of SNPs between all pairs of SNPs
pairwise_nSNP_dist <- pairwise_LD_snp_info %>%
  map(~mutate(., marker2 = factor(marker2, levels = unique(marker2)), 
              marker1 = factor(marker1, levels = levels(marker2)), 
              marker1 = as.numeric(marker1), 
              marker2 = as.numeric(marker2), 
              distance_nSNP = abs(marker1 - marker2)))

# Next calculate the average LD between SNPs that are x = 1, 2, ..., N SNPs apart
pairwise_nSNP_dist_LD <- pairwise_nSNP_dist %>% 
  map_df(~filter(., distance_nSNP <= 1000)) %>% 
  mutate(chrom = as.factor(chrom)) %>% 
  group_by(chrom, distance_nSNP) %>% 
  summarize(meanLD = mean(LD))

# Save the results
save_file <- file.path(result_dir, "S2TP_marker_LD_analysis.RData")
save("pairwise_bpDist_LD", "pairwise_nSNP_dist", "pairwise_LD_snp_info", file = save_file)
