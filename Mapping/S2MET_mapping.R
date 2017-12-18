## Perform GWAS using the S2MET data
## 
## GWAS will be performed using the 1) phenotype BLUEs with environments as a fixed
## effect (i.e. phenotypic mean mapping), 2) all phenotypes to detect QTLxE interaction,


# List of packages to load
# List of packages
packages <- c("dplyr", "purrr", "tibble", "tidyr", "readr", "stringr", "readxl", "parallel", "rrBLUP", "EMMREML")

# packages <- c(packages, "pbr")


# Set the directory of the R packages
package_dir <- NULL
#package_dir <- "/panfs/roc/groups/6/smithkp/neyha001/R/x86_64-pc-linux-gnu-library/3.4/"
package_dir <- "/panfs/roc/groups/6/smithkp/neyha001/R/x86_64-unknown-linux-gnu-library/3.2/"

## Source the function from github
# source("https://raw.github.com/neyhartj/pbr/master/R/gwas.R")
# source("https://raw.github.com/neyhartj/pbr/master/R/gwas_support.R")
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

genos_use <- S2TP_imputed_multi_genos_hmp %>%
  select(rs, chrom, pos, tp_geno)

# Filter the BLUEs to use
S2_MET_BLUEs_use <- S2_MET_BLUEs %>% 
  # filter(!grepl(pattern = "BZI|HTM", x = environment),
  filter(line_name %in% c(tp_geno)) %>%
  group_by(trait, environment) %>%
  filter(sum(line_name %in% tp) > 1) %>%
  ungroup()


# Format the BLUEs across environments
phenos_use <- S2_MET_BLUEs_use %>% 
  select(line_name, environment, trait, value) %>% 
  spread(trait, value)


# Detect cores
n_cores <- detectCores()

## First conduct GWAS of main effects (no QxE) by using the "K" and "G" models
#models <- c("K", "G", "QK", "QG")

#gwas_main <- models %>%
#  map(~gwas(pheno = phenos_use, geno = genos_use, fixed = ~ environment,
#            model = ., impute.method = "pass", n.PC = 2, n.core = n_cores, 
#            test.qxe = FALSE, P3D = TRUE)) %>%
#  set_names(models)

## Save
#save_file <- file.path(result_dir, "S2MET_gwas_genotype_mean.RData")
#save("gwas_main", file = save_file)



## Now run GWAS for the main effect and QxE using a subset of models

# Vector of model types
#models <- c("K", "G", "KE", "GE")
models <- c("G", "GE")


### GWAS of main effect and QxE
# TP only
gwas_main_qxe <- models %>%
  map(~gwas(pheno = phenos_use, geno = genos_use, fixed = ~ environment, 
            model = ., impute.method = "pass", n.PC = 2, n.core = n_cores, 
            test.qxe = TRUE, P3D = TRUE)) %>% 
  set_names(models)

 ## Save
save_file <- file.path(result_dir, "S2MET_gwas_mean_qxe.RData")
save("gwas_main_qxe", file = save_file)


# 
# ## Function for making a matrix full rank
# make_full <- function(X) {
#   svd.X <- svd(X)
#   r <- max(which(svd.X$d > 1e-08))
#   X.full <- as.matrix(svd.X$u[, 1:r])
#   # X.full <- X[,-ncol(X)]
#   return(X.full)
# }
# 
# 
# ### Test out new model that tests QxE using a LRT approach
# pheno_use_test <- phenos_use %>% select(line_name, environment, GrainYield)
# 
# # Marker matrix
# M <- genos_use %>% select(-chrom, -pos) %>% as.data.frame() %>% remove_rownames() %>% column_to_rownames("rs") %>% t()
# 
# # First marker
# snp <- M[,1, drop = FALSE]
# 
# 
# # model.frame
# mf <- model.frame(GrainYield ~ line_name + environment, data = pheno_use_test)
# 
# # y vector
# y <- model.response(mf)
# 
# # Incidence matrix of environments (for subsetting)
# env_inc <- sparse.model.matrix(~ -1 + environment, mf)
# # Vector of mu incidence matrix
# mu <- model.matrix(~ 1, mf)
# 
# # Bind and make full
# X <- cbind(mu, env_inc)
# X_full <- make_full(X = X)
# 
# # Incidence matrix of random genotype effects
# Z <- model.matrix(~ -1 + line_name, mf)
# K1 <- A.mat(X = M, min.MAF = 0, max.missing = 1)
# 
# # Use the env matrix and the snp vector to make a incidence matrix of qxe random effects
# W <- c(Z %*% snp) * env_inc
# K2 <- Diagonal(ncol(W))
# 
# # Add snp main effect to model matrix
# X_use <- cbind(X_full, (Z %*% snp))
#   
# # List of matrices
# rand_list_full <- list(g = list(Z = Z, K = K1),
#                   qxe = list(Z = W, K = K2))
# 
# # For reduced model
# rand_list_red <- list(g = list(Z = Z, K = K1))
# 
# # Fit
# fit_full <- sommer::mmer(Y = y, X = X_use, Z = rand_list_full)
# 
# fit_red <- sommer::mmer(Y = y, X = X_use, Z = rand_list_red)
# 
# # LRT
# lr <- -2 * (fit_red$LL - fit_full$LL)
# 
# 
