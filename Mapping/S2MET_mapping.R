## Perform GWAS using the S2MET data
## 
## GWAS will be performed using the 1) phenotype BLUEs with environments as a fixed
## effect (i.e. phenotypic mean mapping), 2) all phenotypes to detect QTLxE interaction,
## 3) the main effect and stability coefficient from FW regression


# List of packages to load
# List of packages
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

# Extract the tp and vp from the G matrix
s2_imputed_mat_use <- s2_imputed_mat[c(tp_geno, vp_geno),]
# Only the TP
tp_imputed_mat_use <- s2_imputed_mat_use[c(tp_geno),]

# Format the genotype data
s2_imputed_genos_use <- s2_imputed_genos %>% 
  select(marker = `rs#`, chrom, pos, which(names(.) %in% c(tp_geno, vp_geno))) %>%
  as.data.frame()
# Only TP
tp_imputed_genos_use <- s2_imputed_genos_use %>% 
  select(marker:pos, which(names(.) %in% tp_geno))


# Filter the BLUEs to use
S2_MET_BLUEs_use <- S2_MET_BLUEs %>% 
  filter(!grepl(pattern = "BZI|HTM", x = environment),
         line_name %in% c(tp_geno, vp_geno)) %>%
  group_by(trait, environment) %>%
  filter(sum(line_name %in% tp) > 1) %>%
  ungroup()
# Only TP
S2_MET_BLUEs_use_tp <- S2_MET_BLUEs_use %>% 
  filter(line_name %in% tp_geno) %>% 
  droplevels()

# Load the FW results
load(file.path(result_dir, "S2MET_fw_regression_results.RData"))

# Format the phenotypic data
S2_MET_BLUEs_fw_tomodel <- S2_MET_BLUEs_fw %>% 
  filter(line_name %in% c(tp_geno, vp_geno)) %>%
  select(line_name, trait, stability_term, estimate) %>% 
  distinct() %>% 
  spread(stability_term, estimate) %>%
  gather(coef, value, -line_name, -trait) %>% 
  unite("trait_coef", trait, coef) %>% 
  spread(trait_coef, value) %>%
  as.data.frame()

# Format the BLUEs across environments
S2_MET_BLUEs_tomodel <- S2_MET_BLUEs_use %>% 
  select(line_name, environment, trait, value) %>% 
  spread(trait, value)

# Format the BLUEs across environments
S2_MET_BLUEs_tomodel_tp <- S2_MET_BLUEs_use_tp %>% 
  select(line_name, environment, trait, value) %>% 
  spread(trait, value)


# Detect cores
n_cores <- detectCores() 


# Vector of model types
gwas_models <- c("simple", "K", "Q", "QK", "G", "QG")


### GWAS of main effect
# TP only
gwas_tp_main <- map(gwas_models, ~gwas(pheno = S2_MET_BLUEs_tomodel_tp, geno = tp_imputed_genos_use, 
                                       fixed = ~ environment, model = ., impute.method = "pass", 
                                       n.PC = 2, n.core = n_cores))


gwas_tp_main <- set_names(gwas_tp_main, gwas_models)

## Save
save_file <- file.path(result_dir, "S2MET_gwas_tp_model_comparison.RData")
save("gwas_tp_main", file = save_file)


# TP and VP
gwas_all_main <- map(gwas_models, ~gwas(pheno = S2_MET_BLUEs_tomodel, geno = s2_imputed_genos_use, 
                                       fixed = ~ environment, model = ., impute.method = "pass",
                                       n.PC = 2, n.core = n_cores))

gwas_all_main <- set_names(gwas_all_main, gwas_models)

## Save
save_file <- file.path(result_dir, "S2MET_gwas_all_model_comparison.RData")
save("gwas_all_main", file = save_file)





# 
# 
# 
# ## GWAS of genotype main effect and QTLxE
# gwas_qtlxe_out <- gwas_e(pheno = S2_MET_BLUEs_tomodel, geno = s2_imputed_genos_use,
#                          impute.method = "pass", fixed = "environment", n.core = n_cores)
# 
# save_file <- file.path(result_dir, "S2MET_gwas_qtlxe_results.RData")
# save("gwas_qtlxe_out", file = save_file)
# 
# ### GWAS of Genotypic main effect,  and stability across environmental means
# gwas_fw_out <- GWAS(pheno = S2_MET_BLUEs_fw_tomodel, geno = s2_imputed_genos_use, 
#                     n.PC = 0, min.MAF = 0, plot = FALSE, n.core = n_cores)
# 
# 
# ## GWAS of stability coefficients to ECs
# 
# # Load results
# load(file.path(result_dir, "S2MET_ec_fw_regression_results.RData"))
# 
# S2_MET_BLUEs_fw_one_year_tomodel <- S2_MET_BLUEs_one_year_fw %>% 
#   filter(line_name %in% c(tp_geno, vp_geno)) %>% 
#   select(line_name, trait, variable, stability_term, value = estimate) %>%
#   distinct() %>% 
#   unite("trait_coef", trait, variable, stability_term, sep = "_") %>% 
#   spread(trait_coef, value) %>%
#   as.data.frame()
# 
# ### Run GWAS
# gwas_fw_one_year_out <- GWAS(pheno = S2_MET_BLUEs_fw_one_year_tomodel, 
#                              geno = s2_imputed_genos_use, n.PC = 0, min.MAF = 0, 
#                              plot = FALSE, n.core = n_cores)
# 
# 
# 
# 
# # Save
# save_file <- file.path(result_dir, "S2MET_gwas_qtlxe_results.RData")
# save("gwas_qtlxe_out", "gwas_fw_out", "gwas_fw_one_year_out", file = save_file)
