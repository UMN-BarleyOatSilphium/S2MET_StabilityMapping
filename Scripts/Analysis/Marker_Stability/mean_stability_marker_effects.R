## Marker effects for trait mean and trait stability
## 
## 
## This script will calculate marker effects for the trait mean per se and the
## trait stability in order to compare genomic regions that contribute to one or both
## 
## 

# Load packages
packages <- c("dplyr", "purrr", "tibble", "tidyr", "readr", "stringr", "readxl", "modelr", 
              "parallel", "purrrlyr", "rrBLUP", "pbr")

# Set the directory of the R packages
package_dir <- NULL

# Load all packages
invisible(lapply(packages, library, character.only = TRUE, lib.loc = package_dir))

## Directories
proj_dir <- "C:/Users/Jeff/Google Drive/Barley Lab/Projects/S2MET_Mapping//"
alt_proj_dir <- "C:/Users/Jeff/Google Drive/Barley Lab/Projects/S2MET"

# Geno, pheno, and enviro data
geno_dir <-  "C:/Users/Jeff/Google Drive/Barley Lab/Projects/Genomics/Genotypic_Data/GBS_Genotype_Data/"
bopa_geno_dir <- "C:/Users/Jeff/Google Drive/Barley Lab/Projects/Genomics/Genotypic_Data/BOPA_Genotype_Data/"
pheno_dir <- file.path(alt_proj_dir, "Phenotype_Data/")
env_var_dir <- file.path(alt_proj_dir, "Environmental_Variables")

# Other directories
fig_dir <- file.path(proj_dir, "Figures/")
map_dir <- file.path(proj_dir, "Mapping")
entry_dir <- file.path(alt_proj_dir, "Plant_Materials")
analysis_dir <- file.path(proj_dir, "Analysis")
result_dir <- file.path(proj_dir, "Results")

# Load the genotypic data
# load(file.path(geno_dir, "S2_genos_mat.RData"))
# load(file.path(geno_dir, "S2_genos_hmp.RData"))
load(file.path(bopa_geno_dir, "S2TP_multi_genos.RData"))

# Load the FW regression results
load(file.path(result_dir, "S2MET_pheno_mean_fw_results.RData"))

# Load an entry file
entry_list <- read_excel(file.path(entry_dir, "S2MET_project_entries.xlsx"))


# Grab the entry names that are not checks
tp <- entry_list %>% 
  filter(Class == "S2TP") %>% 
  pull(Line)

vp <- entry_list %>% 
  filter(Class == "S2C1R") %>% 
  pull(Line)

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

## ONly use the TP to map

# # Format the genotype data
# genos_use <- s2_imputed_genos %>% 
#   select(marker = `rs#`, chrom, pos, which(names(.) %in% c(tp_geno))) %>%
#   as.data.frame()
genos_use_mat <- S2TP_imputed_multi_genos_mat[tp_geno,]

# Visualize



# Log transform the delta measure
S2_MET_pheno_mean_fw1 <- S2_MET_pheno_mean_fw %>% 
  mutate(delta = log(delta))


# Calculate marker effects for trait mean and trait stability
trait_mean_stab_mar_eff <- S2_MET_pheno_mean_fw1 %>% 
  distinct(trait, line_name, g, b, delta) %>%
  gather(coef, value, -line_name, -trait) %>% # Tidy
  group_by(trait, coef) %>%
  do({
    
    # Grab the data
    df <- .
    
    # Create a model frame and matrices
    mf <- model.frame(value ~ line_name, df)
    y <- model.response(mf)
    Z_use <- model.matrix(~ -1 + line_name, mf)
    
    # Format the marker matrix for the random effect model.matrix
    Z <- Z_use %*% genos_use_mat

    # Fit the model
    fit <- mixed.solve(y = y, Z = Z)
    
    fit$u %>% 
      data.frame(marker = names(.), effect = ., row.names = NULL, stringsAsFactors = FALSE)

  }) %>% ungroup()

# Re-calculate after scaling and centering the values
# Calculate marker effects for trait mean and trait stability
trait_mean_stab_mar_eff_stand <- S2_MET_pheno_mean_fw1 %>% 
  distinct(trait, line_name, g, b, delta) %>%
  gather(coef, value, -line_name, -trait) %>% # Tidy
  group_by(trait, coef) %>%
  mutate(value = scale(value)) %>%
  do({
    
    # Grab the data
    df <- .
    
    # Create a model frame and matrices
    mf <- model.frame(value ~ line_name, df)
    y <- model.response(mf)
    Z_use <- model.matrix(~ -1 + line_name, mf)
    
    # Format the marker matrix for the random effect model.matrix
    Z <- Z_use %*% genos_use_mat
    
    # Fit the model
    fit <- mixed.solve(y = y, Z = Z)
    
    fit$u %>% 
      data.frame(marker = names(.), effect = ., row.names = NULL, stringsAsFactors = FALSE)
    
  }) %>% ungroup()
    
    
# Save this
save_file <- file.path(result_dir, "S2MET_pheno_mean_fw_mar_effect.RData")
save("trait_mean_stab_mar_eff", "trait_mean_stab_mar_eff_stand", file = save_file)
