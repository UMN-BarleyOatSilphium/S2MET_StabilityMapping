## Create marker subsets for dissecting phenotypic stability and mean
##
## This script will create different marker subsets, based on ranking of marker effects,
## plastic/stable markers, and evenly-spaced markers
## 
## Author: Jeff Neyhart
## Last modified: June 26, 2018
## 


# Load the source
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

# Load the stability results
load(file.path(result_dir, "pheno_mean_fw_results.RData"))
# Load the marker effect by environment results
load(file.path(result_dir, "marker_by_env_effects.RData"))
# Load the marker reaction norm results
load(file.path(result_dir, "marker_mean_fw_results.RData"))

# Create the M matrix
M <- S2TP_imputed_multi_genos_mat


## Define the number of markers in the subsets
n_marker_subsets <- c(seq(100, 1000, 100), seq(2000, 9000, 1000)) %>%
  set_names(., str_c("nmar", .))
# Number of subset samples
n_subset_samples <- 100


## First calculate marker effects for the mean, linear stability, and non-linear stability
pheno_mean_fw1 <- pheno_mean_fw %>% 
  distinct(trait, line_name, g, b, delta) %>%
  mutate(log_delta = log(delta)) %>%
  select(-delta)

# Calculate
pheno_mean_fw_marker_effects_all <- pheno_mean_fw1 %>%
  gather(coef, value, g:log_delta) %>% 
  group_by(trait, coef) %>% 
  do(data.frame(marker = colnames(M), effect = mixed.solve(y = .$value, Z = M)$u)) %>%
  ungroup()


## For each marker subset level, drop the markers with the lowest effect, then recalculate effects and repeat
# Reverse the subsets
marker_subsets_sort <- sort(n_marker_subsets, decreasing = TRUE)

# Iterate over the subsets
top_rank_markers <- vector("list", length(marker_subsets_sort))
# Copy the marker effect data.frame
marker_effects <- pheno_mean_fw_marker_effects_all %>%
  mutate(marker = as.character(marker))

for (i in seq_along(marker_subsets_sort)) {
  # Number of markers to keep
  nmar <- marker_subsets_sort[i]
  
  marker_effects_keep <- marker_effects %>% 
    group_by(trait, coef) %>%
    do(arrange(., desc(abs(effect)))[seq(nmar),]) %>% # Ties are resolved arbitrarily
    ungroup() %>%
    left_join(., snp_info, by = "marker") %>%
    select(trait, coef, marker, chrom:cM_pos, effect) %>%
    arrange(trait, coef, chrom, pos)
  
  # Save
  top_rank_markers[[i]] <- marker_effects_keep
  # Add the name
  names(top_rank_markers)[i] <- names(marker_subsets_sort)[i]
  
  # Re-estimate marker effects
  marker_effects <- pheno_mean_fw1 %>%
    gather(coef, value, g:log_delta) %>% 
    group_by(trait, coef) %>% 
    do({
      df <- .
      markers <- subset(marker_effects_keep, trait == unique(df$trait) & coef == unique(df$coef), marker, drop = TRUE)
      data.frame(marker = markers, effect = mixed.solve(y = df$value, Z = M[,markers])$u, row.names = NULL, stringsAsFactors = FALSE)
    }) %>% ungroup()
  
}


# Reverse the order
top_rank_markers <- rev(top_rank_markers)



# Split the markers by chromosome
snp_info_split <- snp_info %>%
  split(.$chrom)

# Number of markers per chromosome
nmar_chr <- map_dbl(snp_info_split, nrow)
# Proportion of markers out of the total
pmar_chr <- nmar_chr / sum(nmar_chr)


## Generate marker intervals for each subset
marker_intervals <- lapply(X = n_marker_subsets, FUN = function(nmar) {
  
  # Weight the number of markers taken by each chromosome by the number of markers on that chromosome
  n_sub_chr <- round(pmar_chr * nmar)
  # Sum
  sum_sub_chr <- sum(n_sub_chr)
  
  # If this is above the total, drop a marker from the chromosome with the most markers
  if (sum_sub_chr > nmar) {
    n_rm <- sum_sub_chr - nmar
    # Which chromosomes to remove
    which_rm <- order(n_sub_chr, decreasing = TRUE)[n_rm]
    # Remove markers
    n_sub_chr[which_rm] <- n_sub_chr[which_rm] - 1
    
  } else if (sum_sub_chr < nmar) {
    # If the sume is less than the total, add a marker to the chromosome with the fewest
    n_add <- nmar - sum_sub_chr
    # Which chromosomes to remove
    which_add <- order(n_sub_chr, decreasing = FALSE)[n_add]
    # Remove markers
    n_sub_chr[which_add] <- n_sub_chr[which_add] + 1
    
  }
  
  list(snp_info_split, n_sub_chr) %>% 
    pmap(~mutate(.x, group = sort(rep(seq(.y), length.out = nrow(.x)))) %>% split(.$group)) 
  
})



## Now use the number of markers in each subset to create evenly-spaced marker groups
# Iterate over the subset number
evenly_spaced_markers <- marker_intervals %>%
  map(., ~map_df(., ~map_df(., ~.[ceiling(mean(c(1, nrow(.)))),])))
    



## For each marker subset level, break the markers into evenly-spaced intervals, take the markers with the lowest effect,
## then re-estimate all marker effects, and repeat

# Iterate over the subsets
top_rank_evenly_spaced_markers <- vector("list", length(marker_subsets_sort))
# Copy the marker effect data.frame
marker_effects <- pheno_mean_fw_marker_effects_all %>%
  mutate(marker = as.character(marker))

# Reverse the marker intervals
marker_intervals_rev <- rev(marker_intervals)

for (i in seq_along(marker_subsets_sort)) {
  # Number of markers to keep
  nmar <- marker_subsets_sort[i]
  
  # Subset the intervals, then collapse the list (chromosomes and groups are retained)
  intervals_use <- marker_intervals_rev[[i]] %>%
    map_df(bind_rows)
  
  # Find the marker in each interval with the largest effect
  marker_effects_keep <- marker_effects %>% 
    group_by(trait, coef) %>%
    do({
      df <- .
      
      # Combine marker effects with the intervals and find the marker with the greatest absolute effect 
      left_join(intervals_use, df, by = "marker") %>% 
        group_by(chrom, group) %>% 
        arrange(desc(abs(effect))) %>% 
        slice(1) %>%
        ungroup()
    }) %>% ungroup() %>%
    select(trait, coef, names(.)) %>%
    arrange(trait, coef, chrom, pos)
  
  # Save
  top_rank_evenly_spaced_markers[[i]] <- marker_effects_keep
  # Add the name
  names(top_rank_evenly_spaced_markers)[i] <- names(marker_subsets_sort)[i]
  
  # Re-estimate marker effects
  marker_effects <- pheno_mean_fw1 %>%
    gather(coef, value, g:log_delta) %>% 
    group_by(trait, coef) %>% 
    do({
      df <- .
      markers <- subset(marker_effects_keep, trait == unique(df$trait) & coef == unique(df$coef), marker, drop = TRUE)
      data.frame(marker = markers, effect = mixed.solve(y = df$value, Z = M[,markers])$u, row.names = NULL, stringsAsFactors = FALSE)
    }) %>% ungroup()
  
}


# Reverse the list
top_rank_evenly_spaced_markers <- rev(top_rank_evenly_spaced_markers)



## Now create random subsets of markers
random_markers <- n_marker_subsets %>%
  map(~replicate(n = n_subset_samples, {
    sample_n(tbl = snp_info, size = .) %>% arrange(chrom, pos) },
    simplify = FALSE))


## Now create markers based on their reaction norms
top_plastic_markers <- marker_mean_fw %>% 
  ungroup() %>% 
  distinct(trait, marker, a, c) %>% 
  gather(coef, effect, a:c) %>%
  left_join(., snp_info) %>%
  group_by(trait, coef) %>%
  do(marker_subsets = {
    df <- .
    
    # Reorder
    df1 <- df %>% arrange(desc(abs(effect)))
    
    n_marker_subsets %>%
      map(~df1[seq(.), ]) %>%
      map(~arrange(., chrom, pos))
    
  }) %>% ungroup()


## Save everything
save("evenly_spaced_markers", "top_rank_evenly_spaced_markers", "top_rank_markers", "top_plastic_markers",
     "random_markers", file = file.path(result_dir, "marker_subsets.RData"))














## Use the training data from TP and VP to create the marker subsets

# Round to the lowest thousand
n_upper <- round(round(ncol(M), digits = -3) - (ncol(M) - round(ncol(M), digits = -3)), -3)

## Define the number of markers in the subsets
n_marker_subsets <- c(seq(100, 1000, 100), seq(2000, n_upper, 1000)) %>%
  set_names(., str_c("nmar", .))
# Number of subset samples
n_subset_samples <- 100

M <- s2_imputed_mat[tp_geno1,]

## First calculate marker effects for the mean, linear stability, and non-linear stability
pheno_mean_fw1 <- pheno_mean_fw_tpvp %>% 
  filter(line_name %in% tp_geno1) %>%
  distinct(trait, line_name, g, b, delta) %>%
  mutate(log_delta = log(delta)) %>%
  select(-delta)

# Calculate
pheno_mean_fw_marker_effects_all <- pheno_mean_fw1 %>%
  gather(coef, value, g:log_delta) %>% 
  group_by(trait, coef) %>% 
  do(data.frame(marker = colnames(M), effect = mixed.solve(y = .$value, Z = M)$u)) %>%
  ungroup()


## For each marker subset level, drop the markers with the lowest effect, then recalculate effects and repeat
# Reverse the subsets
marker_subsets_sort <- sort(n_marker_subsets, decreasing = TRUE)

# Iterate over the subsets
top_rank_markers <- vector("list", length(marker_subsets_sort))
# Copy the marker effect data.frame
marker_effects <- pheno_mean_fw_marker_effects_all %>%
  mutate(marker = as.character(marker))

for (i in seq_along(marker_subsets_sort)) {
  # Number of markers to keep
  nmar <- marker_subsets_sort[i]
  
  marker_effects_keep <- marker_effects %>% 
    group_by(trait, coef) %>%
    do(arrange(., desc(abs(effect)))[seq(nmar),]) %>% # Ties are resolved arbitrarily
    ungroup() %>%
    left_join(., snp_info, by = "marker") %>%
    select(trait, coef, marker, chrom:cM_pos, effect) %>%
    arrange(trait, coef, chrom, pos)
  
  # Save
  top_rank_markers[[i]] <- marker_effects_keep
  # Add the name
  names(top_rank_markers)[i] <- names(marker_subsets_sort)[i]
  
  # Re-estimate marker effects
  marker_effects <- pheno_mean_fw1 %>%
    gather(coef, value, g:log_delta) %>% 
    group_by(trait, coef) %>% 
    do({
      df <- .
      markers <- subset(marker_effects_keep, trait == unique(df$trait) & coef == unique(df$coef), marker, drop = TRUE)
      data.frame(marker = markers, effect = mixed.solve(y = df$value, Z = M[,markers])$u, row.names = NULL, stringsAsFactors = FALSE)
    }) %>% ungroup()
  
}


# Reverse the order
top_rank_markers <- rev(top_rank_markers)



# Split the markers by chromosome
snp_info1_split <- snp_info1 %>%
  split(.$chrom)

# Number of markers per chromosome
nmar_chr <- map_dbl(snp_info1_split, nrow)
# Proportion of markers out of the total
pmar_chr <- nmar_chr / sum(nmar_chr)


## Generate marker intervals for each subset
marker_intervals <- lapply(X = n_marker_subsets, FUN = function(nmar) {
  
  # Weight the number of markers taken by each chromosome by the number of markers on that chromosome
  n_sub_chr <- round(pmar_chr * nmar)
  # Sum
  sum_sub_chr <- sum(n_sub_chr)
  
  # If this is above the total, drop a marker from the chromosome with the most markers
  if (sum_sub_chr > nmar) {
    n_rm <- sum_sub_chr - nmar
    # Which chromosomes to remove
    which_rm <- order(n_sub_chr, decreasing = TRUE)[n_rm]
    # Remove markers
    n_sub_chr[which_rm] <- n_sub_chr[which_rm] - 1
    
  } else if (sum_sub_chr < nmar) {
    # If the sume is less than the total, add a marker to the chromosome with the fewest
    n_add <- nmar - sum_sub_chr
    # Which chromosomes to remove
    which_add <- order(n_sub_chr, decreasing = FALSE)[n_add]
    # Remove markers
    n_sub_chr[which_add] <- n_sub_chr[which_add] + 1
    
  }
  
  list(snp_info1_split, n_sub_chr) %>% 
    pmap(~mutate(.x, group = sort(rep(seq(.y), length.out = nrow(.x)))) %>% split(.$group)) 
  
})



## Now use the number of markers in each subset to create evenly-spaced marker groups
# Iterate over the subset number
evenly_spaced_markers <- marker_intervals %>%
  map(., ~map_df(., ~map_df(., ~.[ceiling(mean(c(1, nrow(.)))),])))




## For each marker subset level, break the markers into evenly-spaced intervals, take the markers with the lowest effect,
## then re-estimate all marker effects, and repeat

# Iterate over the subsets
top_rank_evenly_spaced_markers <- vector("list", length(marker_subsets_sort))
# Copy the marker effect data.frame
marker_effects <- pheno_mean_fw_marker_effects_all %>%
  mutate(marker = as.character(marker))

# Reverse the marker intervals
marker_intervals_rev <- rev(marker_intervals)

for (i in seq_along(marker_subsets_sort)) {
  # Number of markers to keep
  nmar <- marker_subsets_sort[i]
  
  # Subset the intervals, then collapse the list (chromosomes and groups are retained)
  intervals_use <- marker_intervals_rev[[i]] %>%
    map_df(bind_rows)
  
  # Find the marker in each interval with the largest effect
  marker_effects_keep <- marker_effects %>% 
    group_by(trait, coef) %>%
    do({
      df <- .
      
      # Combine marker effects with the intervals and find the marker with the greatest absolute effect 
      left_join(intervals_use, df, by = "marker") %>% 
        group_by(chrom, group) %>% 
        arrange(desc(abs(effect))) %>% 
        slice(1) %>%
        ungroup()
    }) %>% ungroup() %>%
    select(trait, coef, names(.)) %>%
    arrange(trait, coef, chrom, pos)
  
  # Save
  top_rank_evenly_spaced_markers[[i]] <- marker_effects_keep
  # Add the name
  names(top_rank_evenly_spaced_markers)[i] <- names(marker_subsets_sort)[i]
  
  # Re-estimate marker effects
  marker_effects <- pheno_mean_fw1 %>%
    gather(coef, value, g:log_delta) %>% 
    group_by(trait, coef) %>% 
    do({
      df <- .
      markers <- subset(marker_effects_keep, trait == unique(df$trait) & coef == unique(df$coef), marker, drop = TRUE)
      data.frame(marker = markers, effect = mixed.solve(y = df$value, Z = M[,markers])$u, row.names = NULL, stringsAsFactors = FALSE)
    }) %>% ungroup()
  
}


# Reverse the list
top_rank_evenly_spaced_markers <- rev(top_rank_evenly_spaced_markers)



## Now create random subsets of markers
random_markers <- n_marker_subsets %>%
  map(~replicate(n = n_subset_samples, {
    sample_n(tbl = snp_info1, size = .) %>% arrange(chrom, pos) },
    simplify = FALSE))


## Now create markers based on their reaction norms
top_plastic_markers <- marker_mean_fw %>% 
  ungroup() %>% 
  distinct(trait, marker, a, c) %>% 
  gather(coef, effect, a:c) %>%
  left_join(., snp_info) %>%
  group_by(trait, coef) %>%
  do(marker_subsets = {
    df <- .
    
    # Reorder
    df1 <- df %>% arrange(desc(abs(effect)))
    
    n_marker_subsets %>%
      map(~df1[seq(.), ]) %>%
      map(~arrange(., chrom, pos))
    
  }) %>% ungroup()


## Rename everything
evenly_spaced_markers_tpvp <- evenly_spaced_markers
top_rank_evenly_spaced_markers_tpvp <- top_rank_evenly_spaced_markers
top_rank_markers_tpvp <- top_rank_markers
random_markers_tpvp <- random_markers


## Save everything
save("evenly_spaced_markers_tpvp", "top_rank_evenly_spaced_markers_tpvp", "top_rank_markers_tpvp", "random_markers_tpvp",
     file = file.path(result_dir, "marker_subsets_tpvp.RData"))





