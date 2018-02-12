## Find the proportion of GxE variance explained by stable and non-significant markers
## 
## 

# Capture the trait to use (argument 1)
args <- commandArgs(trailingOnly = T)
tr <- args[1]

# List of packages to load
# List of packages
packages <- c("dplyr", "purrr", "tibble", "tidyr", "readr", "stringr", "readxl", "parallel", "rrBLUP",
              "purrrlyr", "sommer")

# Set the directory of the R packages
package_dir <- NULL
package_dir <- "/panfs/roc/groups/6/smithkp/neyha001/R/x86_64-pc-linux-gnu-library/3.4/"
package_dir <- "/panfs/roc/groups/6/smithkp/neyha001/R/x86_64-unknown-linux-gnu-library/3.2/"

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


# Load phenotypic data
load(file.path(pheno_dir, "S2_MET_BLUEs.RData"))
# Load the genotypic data
load(file.path(bopa_geno_dir, "S2TP_multi_genos.RData"))
# Load the marker effect stability results
load(file.path(result_dir, "S2MET_marker_mean_fw_results.RData"))

# Load an entry file
entry_list <- read_excel(file.path(entry_dir, "S2MET_project_entries.xlsx"))


# Grab the TP line names
tp <- entry_list %>% 
  filter(Class == "S2TP") %>% 
  pull(Line)

# Find the intersection of TP lines and those that have been genotyped
tp_geno <- intersect(tp, row.names(S2TP_imputed_multi_genos_mat))


# Create a genotype marker matrix
M <- S2TP_imputed_multi_genos_mat[tp_geno,]

S2_MET_BLUEs_use <- S2_MET_BLUEs %>%
  filter(line_name %in% c(tp_geno)) %>%
  mutate(line_name = as.factor(line_name),
         environment = as.factor(environment))

## First group markers into those that are highly stable, highly sensitive, or neither

S2_MET_marker_mean_fw_tidy <- S2_MET_marker_mean_fw %>% 
  distinct(marker, chrom, pos, trait, b, delta) %>%
  mutate(log_delta = log(delta), b = b + 1) %>%
  dplyr::select(-delta) %>% 
  gather(coef, estimate, -marker:-trait)

# What should be the significance level
alpha <- 0.05

# For each trait, calculate empirical thresholds for significance
S2_MET_marker_fw_sig <- S2_MET_marker_mean_fw_tidy %>%
  filter(coef == "b") %>% 
  group_by(trait) %>% 
  # mutate(estimate = scale(estimate)) %>%
  mutate(lower_perc = quantile(estimate, alpha / 2), 
         upper_perc = quantile(estimate, 1 - (alpha / 2))) %>%
  ungroup() %>%
  mutate(significance = case_when(estimate >= upper_perc ~ "sensitive",
                                  estimate <= lower_perc ~ "stable",
                                  TRUE ~ "average"),
         marker_type = if_else(str_detect(marker, "^S"), "GBS", "BOPA"))


## Fit the mixed models to determine variance components

# Number of model fittings
n_iter <- 25
# Detect cores
n_cores <- detectCores()


## Calculate relationship matrices
# Main effect relationship matrix
K <- A.mat(X = M, min.MAF = 0, max.missing = 1)


## Create marker samples and K matrices
# Extract stable markers
stab_mar <- S2_MET_marker_fw_sig %>% 
  filter(trait == tr, significance == "stable") %>%
  pull(marker)
# Extract sensitive  markers
sens_mar <- S2_MET_marker_fw_sig %>% 
  filter(trait == tr, significance == "sensitive") %>%
  pull(marker)

# What is the minimum number of markers to sample?
min_marker <- min(c(length(stab_mar), length(sens_mar)))

# Extract non-sig markers and take samples
non_sig_mar <- S2_MET_marker_fw_sig %>% 
  filter(trait == tr, significance == "average") %>%
  pull(marker)

non_sig_mar_samples <- replicate(n = n_iter, sample(x = non_sig_mar, size = min_marker),
                                 simplify = FALSE)

## Subset marker matrices and create the static relationship matrices
K_stab <- M[,stab_mar,drop = FALSE] %>% A.mat(X = ., min.MAF = 0, max.missing = 1)
K_sens <- M[,sens_mar,drop = FALSE] %>% A.mat(X = ., min.MAF = 0, max.missing = 1)

# List of non-sig K matrices
K_ns_sample <- non_sig_mar_samples %>%
  map(~M[,.,drop = FALSE]) %>%
  map(~A.mat(X = ., min.MAF = 0, max.missing = 1))




# Extract the data to use
pheno_df <- S2_MET_BLUEs_use %>%
  filter(trait == tr) %>%
  # filter(environment %in% sample(unique(.$environment), 3)) %>%
  droplevels() %>%
  mutate_at(vars(line_name, environment), as.factor)
  

## Model frame/matrices
mf <- model.frame(value ~ line_name + environment + std_error, pheno_df) %>%
  mutate(scale_value = as.numeric(scale(value)),
         avg = interaction(line_name, environment), stab = avg, sens = avg)

# Response
y <- mf$value

# Grand mean as fixed effect
X <- model.matrix(~ 1, mf)

# Weights
wts <- mf$std_error^2
# R matrix
R <- Diagonal(x = wts)
  
## Create the random genotype matrix
# Random effect of genotypic main effect
Z_g <- sparse.model.matrix(~ -1 + line_name, mf)
# Random effect of environment
Z_e <- sparse.model.matrix(~ -1 + environment, mf)
# Random effect of GxE
Z_ge <- sparse.model.matrix(~ -1 + avg, droplevels(mf))

K_e <- Diagonal(ncol(Z_e))

## Compute useable K matrices
# Ze tcrossprod
ZeZe <- tcrossprod(Z_e)

# Stable markers
K_stab_use <- tcrossprod((Z_g %*% K_stab), Z_g)
K_stab_use <- K_stab_use * ZeZe
dimnames(x = K_stab_use) <- list(mf$avg, mf$avg)

# rm(K_stab)

# Plastic markers
K_sens_use <- tcrossprod((Z_g %*% K_sens), Z_g)
K_sens_use <- K_sens_use * ZeZe
dimnames(x = K_sens_use) <- list(mf$avg, mf$avg)

# rm(K_sens)

# Sample of random markers
K_ns_use_sample <- K_ns_sample %>%
  map(~tcrossprod((Z_g %*% .), Z_g))
K_ns_use_sample <- K_ns_use_sample %>%
  map(~. * ZeZe)
K_ns_use_sample <- K_ns_use_sample %>%
  map(~`dimnames<-`(., list(mf$avg, mf$avg)))


# Empty list of variance components
var_comp_out <- vector("list", n_iter)


# Create a list in which the random effects go
Z <- list(
  line_name = list(Z = Z_g, K = K),
  environment = list(Z = Z_e, K = K_e),
  gxe_stable = list(Z = Z_ge, K = K_stab_use),
  gxe_plastic = list(Z = Z_ge, K = K_sens_use),
  gxe_avg = list(Z = Z_ge, K = NULL)
)



## Iterate over the samples
for (i in seq_along(K_ns_use_sample)) {
  
  # Extract the K_matrix
  K_ns_use <- K_ns_use_sample[[i]]
  
  # Edit the Z list
  Z$gxe_avg$K <- K_ns_use
  
  
  # Use sommer to fit the model
  fit <- mmer(Y = y, X = X, Z = Z, R = list(res = R))
  
  # Add the variance components to the list
  var_comp_out[[i]] <- fit$var.comp
  
}


# ## Use BGLR
# K_ns_use <- K_ns_split$`(0.976,7]`[[1]]
# 
# # Create the regression function
# ETA <- list(list(X = Z_e, model = "FIXED"),
#             list(K = kronecker(K_e, K), model = "RKHS"),
#             # list(K = diag(ncol(Z_e)), model = "RKHS"),
#             list(K = K_ns_use, model = "RKHS"))
# 
# fit <- BGLR(y = y, ETA = ETA)
# 
# 
# ## Fit the base model first (accounting for G and E only)
# form_base <- value ~ 1 + (1|line_name) + (1|environment)
# fit_base <- relmatLmer(formula = form_base, data = mf, weights = wts,
#                        relmat = list(line_name = K))
# 
# # Add the residuals to the model frame
# mf1 <- mf %>% 
#   mutate(y_star = resid(fit_base))
# 
# 
# ## Split up these samples by cores
# K_ns_split <- K_ns_use_sample %>%
#   split(., cut(seq_along(.), breaks = n_cores))
# 
# 
# 
# # Apply over cores
# var_comp_out <- mclapply(X = K_ns_split, FUN = function(K_ns_core_list) {
#   
#   # Iterate over the K_ns list of matrices
#   map(K_ns_use_sample, function(K_avg_use) {
#     
#     # # Define the formula
#     # form <- value ~ 1 + (1|line_name) + (1|environment) + (1|stab) + (1|sens) + (1|avg) 
#     # 
#     # # Fit the model
#     # fit <- relmatLmer(formula = form, data = mf, weights = wts,
#     #                   relmat = list(line_name = K, avg = K_avg_use, 
#     #                                 stab = K_stab_use, sens = K_sens_use))
#     
#     # Fit the second model
#     form_ext <- y_star ~ 1 + (1|stab) + (1|sens) + (1|avg)
#     fit_ext <- relmatLmer(formula = form_ext, data = mf1, 
#                           relmat = list(avg = K_avg_use, stab = K_stab_use, sens = K_sens_use))
#   
#     # Return the variance components
#     return(VarProp(fit)) }) })
  

# Save the results
save_file <- file.path(result_dir, str_c("S2MET_pheno_mar_fw_varcomp_", tr, ".RData"))
save("var_comp_out", file = save_file)
