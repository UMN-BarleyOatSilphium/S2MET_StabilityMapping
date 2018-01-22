## S2MET GWAS of Stability Coefficients
## 
## Perform association analysis of the FW stability coefficients for regular phenotype
## FW, and for FW on environmental covariables
## 
## 


# List of packages to load
packages <- c("dplyr", "purrr", "tibble", "tidyr", "readr", "stringr", "readxl", "modelr", 
              "parallel", "purrrlyr", "rrBLUP", "EMMREML", "qvalue")

packages <- c(packages, "pbr")

# Set the directory of the R packages
package_dir <- NULL
package_dir <- "/panfs/roc/groups/6/smithkp/neyha001/R/x86_64-pc-linux-gnu-library/3.4/"
package_dir <- "/panfs/roc/groups/6/smithkp/neyha001/R/x86_64-unknown-linux-gnu-library/3.2/"

source("/panfs/roc/home/smithkp/neyha001/R/my_packages/pbr/R/gwas.R")
source("/panfs/roc/home/smithkp/neyha001/R/my_packages/pbr/R/gwas_support.R")



# Load all packages
invisible(lapply(packages, library, character.only = TRUE, lib.loc = package_dir))


## Directories
proj_dir <- "C:/Users/Jeff/Google Drive/Barley Lab/Projects/S2MET_Mapping//"
proj_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET_Mapping/"

alt_proj_dir <- "C:/Users/Jeff/Google Drive/Barley Lab/Projects/S2MET"
alt_proj_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET/"

# Geno, pheno, and enviro data
geno_dir <-  "C:/Users/Jeff/Google Drive/Barley Lab/Projects/Genomics/Genotypic_Data/GBS_Genotype_Data/"
bopa_geno_dir <- "C:/Users/Jeff/Google Drive/Barley Lab/Projects/Genomics/Genotypic_Data/BOPA_Genotype_Data/"
pheno_dir <- file.path(alt_proj_dir, "Phenotype_Data/")
env_var_dir <- file.path(alt_proj_dir, "Environmental_Variables")

geno_dir <- bopa_geno_dir <-  "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/Data/Genos"
pheno_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/Data/Phenos"
env_var_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/Data/Environmental_Data"


# Other directories
fig_dir <- file.path(proj_dir, "Figures/")
map_dir <- file.path(proj_dir, "Mapping")
entry_dir <- file.path(alt_proj_dir, "Plant_Materials")
analysis_dir <- file.path(proj_dir, "Analysis")
result_dir <- file.path(proj_dir, "Results")


# Load the phenotypic data
load(file.path(pheno_dir, "S2_MET_BLUEs.RData"))
# Load the genotypic data
load(file.path(bopa_geno_dir, "S2TP_multi_genos.RData"))

# Load an entry file
entry_list <- read_excel(file.path(entry_dir, "S2MET_project_entries.xlsx"))


# Grab the entry names that are not checks
tp <- entry_list %>% 
  filter(Class == "S2TP") %>% 
  pull(Line)


# Find the tp that are genotyped
tp_geno <- intersect(tp, row.names(S2TP_imputed_multi_genos_mat))

# Define the checks
checks <- entry_list %>% 
  filter(Class == "Check") %>% 
  pull(Line)

entries <- entry_list %>% 
  pull(Line)


# # Filter the genotype data and create a data.frame for use in GWAS
# # Note 'select' can take a character vector as an argument.
genos_use <- S2TP_imputed_multi_genos_hmp %>%
  select(rs, chrom, pos, tp_geno)


# Load the FW results
load(file.path(result_dir, "S2MET_pheno_mean_fw_results.RData"))
# Load the FW sampling results
load(file.path(result_dir, "S2_MET_pheno_fw_resampling.RData"))

## Format the FW results as a phenotype data.frame
## Filter the lines that only have genotype data
## Log-transform the `delta` stability estimates
## Include BLUEs of the mean
pheno_use <- S2_MET_pheno_mean_fw %>% 
  distinct(line_name, trait, g, b, delta) %>% 
  mutate(log_delta = log(delta)) %>%
  select(-delta) %>%
  gather(coef, value, g, b, log_delta) %>% 
  unite("term", trait, coef, sep = "_") %>% 
  spread(term, value)

  
### GWAS
# Detect cores
n_cores <- detectCores()

# Models to use
models <- c("K", "G", "QK", "QG")

## GWAS for pheno FW
gwas_pheno_mean_fw <- vector("list", length(models)) %>%
  set_names(models)

for (model in models) {
  gwas_pheno_mean_fw[[model]] <- gwas(pheno = pheno_use, geno = genos_use, model = model,
                                 n.PC = 2, P3D = TRUE, n.core = n_cores)

}

save_file <- file.path(result_dir, "S2MET_pheno_fw_mean_gwas_results.RData")
save("gwas_pheno_mean_fw", file = save_file)


## Use the FW resampling data to determine the robustness of the mapping results
# First convert to usable phenotype data.frames and output a list
resample_phenos_use <- S2_MET_pheno_sample_fw %>%
  group_by(p, iter) %>%
  do(pheno_use = {
    select(., line_name, trait, b, delta) %>% 
      mutate(log_delta = log(delta)) %>%
      select(-delta) %>% gather(coef, value, b, log_delta) %>% 
      unite("term", trait, coef, sep = "_") %>% 
      spread(term, value) })

# Iterate over the list of data.frames and perform the association
# Then filter for the significant markers

# First assign cores
resample_phenos_use_list <- resample_phenos_use %>% 
  ungroup() %>% 
  mutate(core = sort(rep(seq(n_cores), length.out = nrow(.)))) %>% 
  split(.$core)

# Iterate over the cores in parallel
resample_gwas_sig_out <- mclapply(X = resample_phenos_use_list, FUN = function(core_df) {
  
  # Create a list
  sig_gwas_list <- vector("list", length = nrow(core_df))
  
  for (i in seq(nrow(core_df))) {
    gwas_out <- gwas(pheno = core_df$pheno_use[[i]], geno = genos_use, model = "G",
                     n.PC = 2, P3D = TRUE, n.core = 1)
    
    # Find the significant hits and return
    sig_gwas_list[[i]] <- gwas_out$scores %>% 
      tbl_df %>% 
      filter(term == "main_effect") %>% 
      group_by(trait) %>% 
      mutate(q_value = qvalue(p = p_value)$qvalue) %>% 
      filter(q_value >= 1)
    
  }
  
  # Add the results to the core_df and return
  core_df %>% 
    mutate(sig_out = sig_gwas_list)
  
}, mc.cores = n_cores)
    

# Save the results
save_file <- file.path(result_dir, "S2MET_pheno_fw_gwas_resample_results.RData")
save("resample_gwas_sig_out", file = save_file)




# # Number of permutations
# n_perm <- 250
# 
# ## Use permutation testing to estimate the sampling distribution
# pheno_fw <- S2_MET_pheno_fw %>%
#   select(line_name, trait, stability_term, estimate) %>%
#   distinct() %>%
#   mutate(estimate = if_else(stability_term == "delta", log(estimate), estimate)) %>%
#   filter(line_name %in% tp_geno) %>%
#   unite("trait_term", trait, stability_term, sep = "_")
# 
# 
# # Generate permutations by each trait-stability level
# pheno_fw_perms <- pheno_fw %>% 
#   split(.$trait_term) %>% 
#   map(~permute(data = ., n = n_perm, line_name)) %>%
#   map("perm") %>%
#   transpose() %>%
#   map(~map_df(., as_data_frame))
# 
# 
# # Re-format for use in the gwas function
# pheno_fw_use <- pheno_fw_perms %>%
#   map(~spread(., trait_term, estimate))
# 
# 
# ## Parallel apply over the list of permutation data.frames
# # First split the list by cores
# pheno_fw_split <- split(pheno_fw_use, cut(x = seq_along(pheno_fw_use), breaks = n_cores))
# 
# 
# # Iterate over the list and parallelize
# pheno_fw_permute_out <- mclapply(X = pheno_fw_split, FUN = function(pheno_fw_perm_list) {
#   
#   # Map over the list of data.frames and run GWAS
#   gwas_out <- pheno_fw_perm_list %>%
#     map(~gwas(pheno = ., geno = genos_use, model = "G", n.PC = 0, P3D = TRUE, n.core = 1))
#   
#   # Extract scores and return a list
#   gwas_out %>% 
#     map("scores") %>% 
#     map(~filter(., term == "main_effect") %>% unnest(estimate))
#   
# }, mc.cores = n_cores)
# 
# 
# # Save the results
# save_file <- file.path(result_dir, str_c("S2MET_pheno_fw_gwas_permutation", n_perm, ".RData"))
# save("pheno_fw_permute_out", file = save_file)
