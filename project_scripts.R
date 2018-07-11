## S2MET GWAS script organizer
## 
## This script includes information on the order in which scripts in this project should be executed.
## It also will check the datestamp of RData files to see check what analyses have been run recently
##

# Repository directory
repo_dir <- getwd()
# Project and other directories
source(file.path(repo_dir, "source.R"))


## Order of scripts
## Below, R code is provided that will run scripts in their entirety. The order of the scripts is the
## recommended order for re-running the analysis

## Part 1 - Phenotypic data analysis
# source(file.path(script_dir, "phenotype_data_summary.R"))   # Summarizes phenotype data across environments

## Part 2 - Calculating genotype mean and stability
# source(file.path(script_dir, "PhenotypicStability/trait_mean_stability_computation.R"))  # Calculates trait mean and stability estimates
# source(file.path(script_dir, "PhenotypicStability/marker_effect_by_env_and_stability.R"))    # Calculates marker effect stability
# source(file.path(script_dir, "PhenotypicStability/marker_subset_generation.R"))    # Generates different marker subsets for predictions (may be removed)

## Part 3 - Association Mapping
# source(file.path(script_dir, "AssociationMapping/stability_mapping.R"))   # Main script for association mapping and MLMM
# source(file.path(script_dir, "AssociationMapping/stability_mapping_resampling.R"))    # Association mapping using environmental subsets (required HPC resources)
# source(file.path(script_dir, "AssociationMapping/pleiotropy_mapping.R"))    # Pleiotropy association mapping






## Check dates of RData files
rdata_files <- list.files(result_dir, pattern = ".RData", full.names = TRUE)
# Load the last file log
load("results_file_log.RData")

# Extract the file information
file_log_new <- rdata_files %>% 
  map_df(~data_frame(file = basename(.), info = list(file.info(.)))) %>% 
  unnest() %>% 
  select(file, change_time = ctime)

# Compare against the previous file log and output the files that have been changed
file_log %>% 
  left_join(., rename(file_log_new, change_time_new = change_time), by = "file") %>% 
  filter(change_time_new > change_time)





