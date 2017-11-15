## S2MET GWAS of Stability Coefficients
## 
## Perform association analysis of the FW stability coefficients for regular phenotype
## FW, and for FW on environmental covariables
## 
## 


# List of packages to load
packages <- c("dplyr", "purrr", "tibble", "tidyr", "readr", "stringr", "readxl", "modelr", 
              "parallel", "purrrlyr", "rrBLUP", "pbr")

# Set the directory of the R packages
package_dir <- NULL
package_dir <- "/panfs/roc/groups/6/smithkp/neyha001/R/x86_64-pc-linux-gnu-library/3.4/"

# Load all packages
invisible(lapply(packages, library, character.only = TRUE, lib.loc = package_dir))


## Directories
proj_dir <- "C:/Users/Jeff/Google Drive/Barley Lab/Projects/S2MET_Mapping//"
proj_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET_Mapping/" 

alt_proj_dir <- "C:/Users/Jeff/Google Drive/Barley Lab/Projects/S2MET"
alt_proj_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET/"

# Geno, pheno, and enviro data
geno_dir <-  "C:/Users/Jeff/Google Drive/Barley Lab/Projects/Genomics/Genotypic_Data/GBS_Genotype_Data/"
pheno_dir <- file.path(alt_proj_dir, "Phenotype_Data/")
env_var_dir <- file.path(alt_proj_dir, "Environmental_Variables")

geno_dir <-  "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/Data/GBS_Genos"
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
load(file.path(geno_dir, "S2_genos_mat.RData"))
load(file.path(geno_dir, "S2_genos_hmp.RData"))
# Load environmental data
load(file.path(env_var_dir, "environmental_data_compiled.RData"))

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
tp_geno <- intersect(tp, row.names(s2_imputed_mat))
vp_geno <- intersect(vp, row.names(s2_imputed_mat))

# Define the checks
checks <- entry_list %>% 
  filter(Class == "Check") %>% 
  pull(Line)

entries <- entry_list %>% 
  pull(Line)


## List of populations to use
pops <- list(pop_all = c(tp_geno, vp_geno), pop_tp = tp_geno)

# Filter the genotype data and create data.frames for use in GWAS
# Note 'select' can take a character vector as an argument.
genos_use <- pops %>%
  map(~select(s2_imputed_genos, marker = `rs#`, chrom, pos, .))
  

# Load the FW results
load(file.path(result_dir, "S2MET_pheno_fw_regression_results.RData"))
load(file.path(result_dir, "S2MET_ec_fw_regression_results.RData"))

## Format the FW results as a phenotype data.frame
## Filter the lines that only have genotype data
pheno_fw_use <- list(S2_MET_pheno_fw, pops) %>%
  pmap(~select(.x, line_name, trait, stability_term, estimate) %>%
        distinct() %>%
        filter(line_name %in% .y) %>%
        unite("trait_term", trait, stability_term, sep = "_") %>% 
        spread(trait_term, estimate) )

ec_oneyear_fw_use <- list(S2_MET_ec_oneyear_fw, pops) %>% 
  pmap(~select(.x, line_name, trait, variable, stability_term, estimate) %>% 
        distinct() %>% 
        filter(line_name %in% .y) %>%
        unite(trait_variable_term, trait, variable, stability_term, sep = "_") %>% 
        spread(trait_variable_term, estimate) )

ec_multiyear_fw_use <- list(S2_MET_ec_multiyear_fw, pops) %>% 
  pmap(~select(.x, line_name, trait, variable, stability_term, estimate) %>% 
         distinct() %>% 
         filter(line_name %in% .y) %>%
         unite(trait_variable_term, trait, variable, stability_term, sep = "_") %>% 
         spread(trait_variable_term, estimate) )


### GWAS
# Detect cores
n_cores <- detectCores()

# Models to use
models <- c("K", "G")

## GWAS for pheno FW
gwas_pheno_fw <- models %>%
  map(function(mod) 
    list(pheno_fw_use, genos_use) %>% 
      pmap(~gwas(pheno = .x, geno = .y, model = mod, P3D = TRUE, n.core = n_cores)))


save_file <- file.path(result_dir, "S2MET_pheno_fw_gwas_results.RData")
save("gwas_pheno_fw", file = save_file)



# GWAS for single-year ECs
gwas_singleyear_ec_fw <- models %>%
  map(function(mod) 
    list(ec_oneyear_fw_use, genos_use) %>% 
      pmap(~gwas(pheno = .x, geno = .y, model = mod, P3D = TRUE, n.core = n_cores)))

# GWAS for multi-year ECs
gwas_multiyear_ec_fw <- models %>%
  map(function(mod) 
    list(ec_multiyear_fw_use, genos_use) %>% 
      pmap(~gwas(pheno = .x, geno = .y, model = mod, P3D = TRUE, n.core = n_cores)))


## Save
save_file <- file.path(result_dir, "S2MET_fw_gwas_results.RData")
save("gwas_singleyear_ec_fw", "gwas_multiyear_ec_fw", file = save_file)











