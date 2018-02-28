## S2MET Mapping
## Calculate marker effects in each environment
## 
## This script will calculate marker effects across environment using a G model, similar to GWAS.
## 

# # Capture the trait to use (argument 1)
# args <- commandArgs(trailingOnly = T)
# tr <- args[1]

# List of packages to load
# List of packages
packages <- c("dplyr", "purrr", "tibble", "tidyr", "readr", "stringr", "readxl", 
              "Matrix", "lme4qtl", "sommer", "rrBLUP")

# packages <- c(packages, "pbr")

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

bopa_geno_dir <- geno_dir <-  "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/Data/Genos"
pheno_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/Data/Phenos"


# Other directories
script_dir <- file.path(proj_dir, "Scripts/")
analysis_dir <- file.path(script_dir, "Analysis")
fig_dir <- file.path(script_dir, "Figures/")
map_dir <- file.path(analysis_dir, "GWAS")
result_dir <- file.path(proj_dir, "Results")


# Load the phenotypic data
load(file.path(pheno_dir, "S2_MET_BLUEs.RData"))
# Load the genotypic data
load(file.path(bopa_geno_dir, "S2TP_multi_genos.RData"))

# Matrix of genotype data for all individuals
# This must remain intact, and then will be subsetted below
M <- S2TP_imputed_multi_genos_mat

# Get the marker names
markers <- colnames(M)


# Remove the environments in which the vp was only observed
S2_MET_BLUEs_use <- S2_MET_BLUEs %>%
  filter(line_name %in% row.names(M)) %>%
  mutate(line_name = as.factor(line_name),
         environment = as.factor(environment))


# ## Use the GWAS G model to estimate the effect of each marker in each environment
# ## This will be the environment-specific marker effect + the mean
# 
# # The code below using the lme4qtl package
# 
# # Subset markers by chromosome
# markers_by_chrom <- S2TP_imputed_multi_genos_hmp %>%
#   select(marker = rs, chrom)
# 
# # All marker names
# markers <- colnames(M)
# 
# # Create relationship matrices per chromosome
# K_chr <- markers_by_chrom %>%
#   split(.$chrom) %>%
#   map(~M[,setdiff(markers, .$marker),drop = FALSE]) %>%
#   map(~A.mat(X = ., min.MAF = 0, max.missing = 1))
# 
# ## One trait at a time
# trait_df <- S2_MET_BLUEs_use %>%
#   filter(trait == tr) %>%
#   droplevels()
# 
# 
# 
# # Model matrices for line_name and environment
# Zg <- model.matrix(~ -1 + line_name, trait_df)
# Ze <- model.matrix(~ -1 + environment, trait_df)
# 
# # Extract weights
# wts <- trait_df$std_error^2
# 
# # Split markers by core
# core_list <- markers_by_chrom %>%
#   mutate(core = sort(rep(seq(n_cores), length.out = nrow(.)))) %>%
#   split(.$core)
# 
# # Iterate over the core list
# marker_score_out <- mclapply(X = core_list, FUN = function(core) {
# 
#   # empty list to store results
#   core_list_out <- vector("list", nrow(core)) %>%
#     set_names(core$marker)
# 
#   # Apply a function over the marker matrix
#   for (i in seq(nrow(core))) {
#     # Subset the snp
#     snp <- core[i,]
# 
#     mar <- Zg %*% M[,snp$marker, drop = FALSE]
# 
#     K_use <- K_chr[[snp$chrom]]
# 
#     # fit the model
#     fit <- relmatLmer(value ~ mar:environment + (1|line_name) + (1|environment), trait_df,
#                       relmat = list(line_name = K_use), weights = wts)
# 
#     # Extract the coefficients
#     core_list_out[[i]] <- tidy(fit) %>%
#       subset(str_detect(term, "mar"), c(term, estimate))
# 
#   }
# 
#   # return the list
#   return(core_list_out)
# 
# }, mc.cores = n_cores)
# 

# # Save the output
# save_file <- file.path(result_dir, str_c("S2MET_marker_by_env_", tr, ".RData"))
# save("marker_score_out", file = save_file)


## Calculate marker effects separately for each environment
marker_by_env <- S2_MET_BLUEs_use %>%
  group_by(trait, environment) %>%
  do({
    df <- .
    
    # Model frame
    mf <- model.frame(value ~ line_name + std_error, df)

    # Response
    y <- model.response(mf)
    # Grand mean design
    X <- model.matrix(~ 1, mf)
    # Line_name design
    Zg <- model.matrix(~ -1 + line_name, mf)
    
    # Subset the marker matrix
    Z <- Zg %*% M
    K <- diag(ncol(Z))
    
    # Pull out the weights into an R matrix
    # R <- solve(diag(mf$std_error^2))
    
    # Fit
    # fit <- sommer::mmer(Y = y, X = X, Z = list(marker = list(Z = Z, K = K)), R = list(res = R))
    fit <- mixed.solve(y = y, Z = Z, method = "REML")
    
    
    # # Grab the marker effects and return
    # marker_effect <- fit$u.hat$marker %>% 
    #   as.data.frame() %>% 
    #   rownames_to_column("marker") %>% 
    #   dplyr::select(marker, effect = T1)
    
    # Grab the marker effects and return
    marker_effect <- fit$u %>% 
      as.data.frame() %>% 
      rownames_to_column("marker") %>%
      rename(effect = ".")
    
    # Return a data_frame
    data_frame(marker_effect = list(marker_effect), fit = list(fit)) })
    

# ## Calculate marker effects x environment as in Lopez-Cruz 2015
# 
# # Find the unique traits
# trts <- unique(S2_MET_BLUEs_use$trait)
# 
# # Create an empty list
# marker_by_env <- list()
# 
# # Iterate over traits
# for (tr in trts) {
#   
#   # Extract the data for that trait
#   df <- S2_MET_BLUEs_use %>%
#     filter(trait == tr)
#     
#   # Create a model frame with lines and environments
#   mf <- model.frame(value ~ line_name + environment + std_error, data = df, drop.unused.levels = TRUE)
#   env_names <- levels(mf$environment)
#   
#   # Response vector
#   y <- model.response(mf)
#   
#   # Model matrix for the fixed effects
#   X <- model.matrix(~ environment, data = mf)
#   
#   # Model matrix of the line_name incidence
#   Z_line <- sparse.model.matrix(~ -1 + line_name, mf)
#   
#   # Model matrix of the marker random effects
#   ## First the main effect of each marker
#   Z0 <- Z_line %*% M
#   K0 <- Diagonal(ncol(Z0))
#   
#   # Random deviation of each marker in each environment
#   Z1 <- mf %>%
#     split(.$environment) %>%
#     map(~ sparse.model.matrix(~ -1 + line_name, .)) %>%
#     map(~ . %*% M)
#   Z1 <- .bdiag(Z1)
#   
#   K1 <- Diagonal(ncol(Z1))
# 
#   # fit <- EMMREML::emmremlMultiKernel(y = y, X = X, Zlist = list(Z0, Z1), Klist = list(K0, K1))
#   fit <- mmer(Y = y, X = X, Z = list(m = list(Z = Z0, K = K0), mxe = list(Z = Z1, K = K1)))
#   
#   # Extract the base marker effects
#   base_eff <- fit$u.hat$m %>%
#     as.data.frame() %>%
#     rownames_to_column("marker") %>%
#     mutate(environment = NA) %>%
#     rename(effect = T1)
#   
#   # Extract the mxe effects
#   mxe_eff <- fit$u.hat$mxe %>% 
#     as.data.frame() %>% 
#     mutate(marker = rep(markers, length(env_names)), 
#            environment = rep(env_names, each = length(markers))) %>% 
#     rename(effect = T1)
#   
#   # combine
#   random_eff <- bind_rows(base_eff, mxe_eff)
#   
#   # Add to the list
#   marker_by_env[[tr]] <- random_eff
#   
# }

## Save this
save_file <- file.path(result_dir, "S2MET_marker_ranef_by_env.RData")
save("marker_by_env", file = save_file)

