## S2MET GWAS of Stability Coefficients
## 
## Perform association analysis of the FW stability coefficients for regular phenotype
## FW, and for FW on environmental covariables
## 
## 


# List of packages to load
packages <- c("dplyr", "purrr", "tibble", "tidyr", "readr", "stringr", "readxl", "modelr", 
              "parallel", "purrrlyr", "rrBLUP", "EMMREML")

#packages <- c(packages, "pbr")

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
# load(file.path(geno_dir, "S2_genos_mat.RData"))
# load(file.path(geno_dir, "S2_genos_hmp.RData"))
load(file.path(bopa_geno_dir, "S2TP_multi_genos.RData"))
# Load environmental data
load(file.path(env_var_dir, "environmental_data_compiled.RData"))

# Load an entry file
entry_list <- read_excel(file.path(entry_dir, "S2MET_project_entries.xlsx"))


# Grab the entry names that are not checks
tp <- entry_list %>% 
  filter(Class == "S2TP") %>% 
  pull(Line)

# vp <- entry_list %>% 
#   filter(Class == "S2C1R") %>% 
#   pull(Line)

# Find the tp and vp that are genotypes
# tp_geno <- intersect(tp, row.names(s2_imputed_mat))
tp_geno <- intersect(tp, row.names(S2TP_imputed_multi_genos_mat))
# vp_geno <- intersect(vp, row.names(s2_imputed_mat))

# Define the checks
checks <- entry_list %>% 
  filter(Class == "Check") %>% 
  pull(Line)

entries <- entry_list %>% 
  pull(Line)


# # Filter the genotype data and create a data.frame for use in GWAS
# # Note 'select' can take a character vector as an argument.
# genos_use <- s2_imputed_genos %>%
#   select(marker = `rs#`, chrom, pos, tp_geno)

genos_use <- S2TP_imputed_multi_genos_hmp %>%
  select(rs, chrom, pos, tp_geno)


# Load the FW results
load(file.path(result_dir, "S2MET_pheno_fw_regression_results.RData"))
load(file.path(result_dir, "S2MET_ec_fw_regression_results.RData"))

## Format the FW results as a phenotype data.frame
## Filter the lines that only have genotype data
## Log-transform the `delta` stability estimates
pheno_fw_use <- S2_MET_pheno_fw %>%
  select(line_name, trait, stability_term, estimate) %>%
  distinct() %>%
  mutate(estimate = if_else(stability_term == "delta", log(estimate), estimate)) %>%
  filter(line_name %in% tp_geno) %>%
  unite("trait_term", trait, stability_term, sep = "_") %>% 
  spread(trait_term, estimate)
  

#ec_oneyear_fw_use <- S2_MET_ec_oneyear_fw %>% 
#  select(line_name, trait, variable, stability_term, estimate) %>% 
#     distinct() %>% 
#     filter(line_name %in% tp_geno) %>%
#     unite(trait_variable_term, trait, variable, stability_term, sep = "_") %>% 
#     spread(trait_variable_term, estimate)

#ec_multiyear_fw_use <- S2_MET_ec_multiyear_fw %>% 
#  select(line_name, trait, variable, stability_term, estimate) %>% 
#  distinct() %>%  
#  filter(line_name %in% tp_geno) %>%
#  unite(trait_variable_term, trait, variable, stability_term, sep = "_") %>% 
#  spread(trait_variable_term, estimate)


### GWAS
# Detect cores
n_cores <- detectCores()

# # Models to use
# models <- c("K", "G", "QK", "QG")
# 
# ## GWAS for pheno FW
# gwas_pheno_fw <- vector("list", length(models)) %>%
#   set_names(models)
# 
# for (model in models) {
#   gwas_pheno_fw[[model]] <- gwas(pheno = pheno_fw_use, geno = genos_use, model = model, 
#                                  n.PC = 2, P3D = TRUE, n.core = n_cores)
#   
# }
# 
# # save_file <- file.path(result_dir, "S2MET_pheno_fw_gwas_results.RData")
# save_file <- file.path(result_dir, "S2MET_pheno_fw_gwas_multi_geno_results.RData")
# save("gwas_pheno_fw", file = save_file)



# Number of permutations
n_perm <- 300

## Use permutation testing to estimate the sampling distribution
pheno_fw <- S2_MET_pheno_fw %>%
  select(line_name, trait, stability_term, estimate) %>%
  distinct() %>%
  mutate(estimate = if_else(stability_term == "delta", log(estimate), estimate)) %>%
  filter(line_name %in% tp_geno) %>%
  unite("trait_term", trait, stability_term, sep = "_")


# Generate permutations by each trait-stability level
pheno_fw_perms <- pheno_fw %>% 
  split(.$trait_term) %>% 
  map(~permute(data = ., n = n_perm, line_name)) %>%
  map("perm") %>%
  transpose() %>%
  map(~map_df(., as_data_frame))


# Re-format for use in the gwas function
pheno_fw_use <- pheno_fw_perms %>%
  map(~spread(., trait_term, estimate))


## Parallel apply over the list of permutation data.frames
# First split the list by cores
pheno_fw_split <- split(pheno_fw_use, cut(x = seq_along(pheno_fw_use), breaks = n_cores))


# Iterate over the list and parallelize
pheno_fw_permute_out <- mclapply(X = pheno_fw_split, FUN = function(pheno_fw_perm_list) {
  
  # Map over the list of data.frames and run GWAS
  gwas_out <- pheno_fw_perm_list %>%
    map(~gwas(pheno = ., geno = genos_use, model = "G", n.PC = 0, P3D = TRUE, n.core = 1))
  
  # Extract scores and return a list
  gwas_out %>% 
    map("scores") %>% 
    map(~filter(., term == "main_effect") %>% unnest(estimate))
  
}, mc.cores = n_cores)


# Save the results
save_file <- file.path(result_dir, str_c("S2MET_pheno_fw_gwas_permutation", n_perm, ".RData")
save("pheno_fw_permute_out", file = save_file)


# ### GWAS for stability across ECs
# # GWAS for single-year ECs
# gwas_singleyear_ec_fw <- models %>%
#   map(., ~gwas(pheno = ec_oneyear_fw_use, geno = genos_use, model = ., P3D = TRUE, n.core = n_cores))
# 
# # GWAS for multi-year ECs
# gwas_multiyear_ec_fw <- models %>%
#   map(., ~gwas(pheno = ec_multiyear_fw_use, geno = genos_use, model = ., P3D = TRUE, n.core = n_cores))
# 
# 
# ## Save
# save_file <- file.path(result_dir, "S2MET_fw_gwas_results.RData")
# save("gwas_singleyear_ec_fw", "gwas_multiyear_ec_fw", file = save_file)
# 
# 
# 
# 
# 
# 
# ### Archived code for testing the TP and the TP + VP
# 
# 
# ## List of populations to use
# pops <- list(pop_all = c(tp_geno, vp_geno), pop_tp = tp_geno)
# 
# # Filter the genotype data and create data.frames for use in GWAS
# # Note 'select' can take a character vector as an argument.
# genos_use <- pops %>%
#   map(~select(s2_imputed_genos, marker = `rs#`, chrom, pos, .))
#   
# 
# # Load the FW results
# load(file.path(result_dir, "S2MET_pheno_fw_regression_results.RData"))
# load(file.path(result_dir, "S2MET_ec_fw_regression_results.RData"))
# 
# ## Format the FW results as a phenotype data.frame
# ## Filter the lines that only have genotype data
# pheno_fw_use <- list(S2_MET_pheno_fw, pops) %>%
#   pmap(~select(.x, line_name, trait, stability_term, estimate) %>%
#         distinct() %>%
#         filter(line_name %in% .y) %>%
#         unite("trait_term", trait, stability_term, sep = "_") %>% 
#         spread(trait_term, estimate) )
# 
# ec_oneyear_fw_use <- list(S2_MET_ec_oneyear_fw, pops) %>% 
#   pmap(~select(.x, line_name, trait, variable, stability_term, estimate) %>% 
#         distinct() %>% 
#         filter(line_name %in% .y) %>%
#         unite(trait_variable_term, trait, variable, stability_term, sep = "_") %>% 
#         spread(trait_variable_term, estimate) )
# 
# ec_multiyear_fw_use <- list(S2_MET_ec_multiyear_fw, pops) %>% 
#   pmap(~select(.x, line_name, trait, variable, stability_term, estimate) %>% 
#          distinct() %>% 
#          filter(line_name %in% .y) %>%
#          unite(trait_variable_term, trait, variable, stability_term, sep = "_") %>% 
#          spread(trait_variable_term, estimate) )
# 
# 
# ### GWAS
# # Detect cores
# n_cores <- detectCores()
# 
# # Models to use
# models <- c("K", "G")
# 
# ## GWAS for pheno FW
# gwas_pheno_fw <- models %>%
#   map(function(mod) 
#     list(pheno_fw_use, genos_use) %>% 
#       pmap(~gwas(pheno = .x, geno = .y, model = mod, P3D = TRUE, n.core = n_cores)))
# 
# 
# save_file <- file.path(result_dir, "S2MET_pheno_fw_gwas_results.RData")
# save("gwas_pheno_fw", file = save_file)
# 
# 
# 
# # GWAS for single-year ECs
# gwas_singleyear_ec_fw <- models %>%
#   map(function(mod) 
#     list(ec_oneyear_fw_use, genos_use) %>% 
#       pmap(~gwas(pheno = .x, geno = .y, model = mod, P3D = TRUE, n.core = n_cores)))
# 
# # GWAS for multi-year ECs
# gwas_multiyear_ec_fw <- models %>%
#   map(function(mod) 
#     list(ec_multiyear_fw_use, genos_use) %>% 
#       pmap(~gwas(pheno = .x, geno = .y, model = mod, P3D = TRUE, n.core = n_cores)))
# 
# 
# ## Save
# save_file <- file.path(result_dir, "S2MET_fw_gwas_results.RData")
# save("gwas_singleyear_ec_fw", "gwas_multiyear_ec_fw", file = save_file)
# 
