## Create marker subsets for dissecting phenotypic stability and mean
##
## This script will create different marker subsets, based on ranking of marker effects,
## plastic/stable markers, and evenly-spaced markers
## 
## Author: Jeff Neyhart
## Last modified: June 23, 2018
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

# Calculate marker effects
pheno_mean_fw_marker_effects <- pheno_mean_fw1 %>%
  gather(coef, value, g:log_delta) %>% 
  group_by(trait, coef) %>% 
  do(data.frame(marker = colnames(M), effect = mixed.solve(y = .$value, Z = M)$u)) %>%
  ungroup()

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
    
## Now use the marker effects to find the marker with the greatest absolute value in each
## interval
top_rank_evenly_spaced_markers <- pheno_mean_fw_marker_effects %>%
  left_join(., snp_info) %>%
  group_by(trait, coef) %>%
  do(marker_subsets = {
    df <- . 
    df1 <- df %>% split(.$chrom) %>% map(~select(., -chrom:-cM_pos))
    
    marker_intervals %>%
      map(~list(., df1) %>% pmap(~{
        lapply(X = .x, left_join, y = .y, by = "marker") %>% 
          map_df(~subset(., abs(effect) == max(abs(effect)))[1,]) }) %>%
          bind_rows() ) # Choose the first marker if more than one 
    })

# Ungroup
top_rank_evenly_spaced_markers <- ungroup(top_rank_evenly_spaced_markers)


## Now simply take the markers with the top n effects for each character
top_rank_markers <- pheno_mean_fw_marker_effects %>%
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


## Now create random subsets of markers
random_markers <- n_marker_subsets %>%
  map(~replicate(n = n_subset_samples, {
    sample_n(tbl = snp_info, size = .) %>% arrange(chrom, pos) },
    simplify = FALSE))


## Now create markers based on their reaction norms
top_plastic_markers <- marker_effect_stab %>% 
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











