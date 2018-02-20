## S2MET GWAS of Stability Coefficients
## 
## Perform association analysis of the FW stability coefficients for regular phenotype
## FW, and for FW on environmental covariables
## 
## 


# List of packages to load
packages <- c("dplyr", "purrr", "tibble", "tidyr", "readr", "stringr", "readxl", "modelr", 
              "parallel", "purrrlyr", "rrBLUP", "EMMREML", "qvalue", "lme4qtl")

# packages <- c(packages, "pbr")

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
script_dir <- file.path(proj_dir, "Scripts/")
analysis_dir <- file.path(script_dir, "Analysis")
fig_dir <- file.path(script_dir, "Figures/")
map_dir <- file.path(analysis_dir, "GWAS")
result_dir <- file.path(proj_dir, "Results")


# Load the phenotypic data
load(file.path(pheno_dir, "S2_MET_BLUEs.RData"))
# Load the genotypic data
load(file.path(bopa_geno_dir, "S2TP_multi_genos.RData"))

# Load an entry file
entry_list <- read_excel(file.path(entry_dir, "S2MET_project_entries.xlsx"))

# Find the tp that are genotyped
tp_geno <- row.names(S2TP_imputed_multi_genos_mat)


# # Filter the genotype data and create a data.frame for use in GWAS
# # Note 'select' can take a character vector as an argument.
genos_use <- S2TP_imputed_multi_genos_hmp %>%
  select(rs, chrom, pos, tp_geno)



# Load the FW results
load(file.path(result_dir, "S2MET_pheno_mean_fw_results.RData"))
# Load the FW sampling results
load(file.path(result_dir, "S2MET_pheno_fw_resampling.RData"))


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


### Mapping using the observed stability estimates ###

# Models to use
models <- c("K", "G", "QK", "QG")

## GWAS for pheno FW
gwas_pheno_mean_fw <- vector("list", length(models)) %>%
  set_names(models)

for (model in models) {
  gwas_pheno_mean_fw[[model]] <- gwas(pheno = pheno_use, geno = genos_use, model = model,
                                 n.PC = 2, P3D = TRUE, n.core = n_cores)

}


### Fit a multi-locus mixed model to identify significant QTL ###

# First correct the p_values
gwas_pheno_mean_fw_adj <- gwas_pheno_mean_fw %>%
  map_df(~mutate(.$scores, model = .$metadata$model)) %>%
  filter(term == "main_effect") %>%
  separate(trait, c("trait", "coef"), sep = "_", extra = "merge") %>%
  group_by(., trait, coef, model) %>%
  mutate(p_adj = p.adjust(p_value, method = "fdr"),
         q_value = qvalue(p = p_value)$qvalue,
         local_fdr = qvalue(p = p_value)$lfdr,
         neg_log10_fdr05 = -log10(0.05),
         neg_log10_fdr10 = -log10(0.10)) %>%
  mutate_at(vars(p_value, p_adj, q_value), funs(neg_log10 = -log10(.))) %>%
  ungroup()

# Cutoff for significance
sig_cutoff <- 0.05

# Select the G model and identify the signficant associations
gwas_adj_sig <- gwas_pheno_mean_fw_adj %>%
  filter(model == "G", q_value <= sig_cutoff)


## Fit MLMM on a chromosome-wise basis
# Identify the SNPs on each chromosome
snps_by_chrom <- S2TP_imputed_multi_genos_hmp %>%
  split(.$chrom) %>%
  map("rs")

# Create relationship matrices for each chromosome
K_c <- snps_by_chrom %>%
  map(~subset(S2TP_imputed_multi_genos_mat, select = setdiff(colnames(S2TP_imputed_multi_genos_mat), .))) %>%
  map(~A.mat(., min.MAF = 0, max.missing = 1))

# Define a function to make a matrix full rank
make_full <- function(X) {
  svd.X <- svd(X)
  r <- max(which(svd.X$d > 1e-08))
  X.full <- as.matrix(svd.X$u[, 1:r])
  # X.full <- X[,-ncol(X)]
  return(X.full)
}

# Define a function that fits the linear mm and returns estimates, std.errors, and pvalue of
# markers
mlmm <- function(y, X, Z, K) {
  # Fit the model
  fit <- mixed.solve(y = y, Z = Z, K = K, X = X, SE = TRUE)
  alpha <- fit$beta[-1]
  se <- fit$beta.SE[-1]
  # Perform hypothesis test
  chisq <- alpha^2 / se^2
  p_value <- pchisq(q = chisq, df = 1, lower.tail = FALSE)
  q_value <- p.adjust(p = p_value, method = "BH")

  # Fit a reduced model (no SNPs)
  X_red <- matrix(1, nrow = length(y))
  fit_red <- mixed.solve(y = y, Z = Z, K = K, X = X_red, SE = TRUE)
  
  # Calculate R-squared
  y_hat_red <- ((X_red %*% fit_red$beta) + (Z %*% fit_red$u))
  y_hat_snp <- (X %*% fit$beta) + (Z %*% fit$u)
  
  # Calculate R^2 and adjusted R^2 for the full model
  n <- length(y)
  p <- ncol(X) - 1

  R_sqr <- cor(y, y_hat_snp)^2
  R_sqr_adj <- R_sqr - ((1 - R_sqr) * (p / (n - p - 1)))
  
  # Now calculate R^2 and adjusted R^2 for the reduced model
  n <- length(y)
  p <- ncol(X_red) - 1
  
  R_sqr_red <- cor(y, y_hat_red)^2
  R_sqr_red_adj <- R_sqr_red - ((1 - R_sqr_red) * (p / (n - p - 1)))

  # Return a list and a data.frame
  data.frame(marker = names(alpha), alpha, se, p_value, q_value, R_sqr = R_sqr,
             R_sqr_adj = R_sqr_adj)

}

# Function for removing colinear SNPs
remove_colinear <- function(x) {
  # LD
  x_LD <- cor(x)^2
  # Find the index of the SNPs to be dropped
  snp_index <- apply(X = x_LD, MARGIN = 1, FUN = function(LD) max(which(LD > 0.99)))
  to_drop <- which(duplicated(snp_index))

  ifelse(length(to_drop) > 0, return(x[,-to_drop, drop = FALSE]), return(x))

}

# Iterate over trait, coef, and chromosome
gwas_mlmm_Gmodel <- gwas_adj_sig %>%
  group_by(trait, coef, chrom) %>%
  do({
    # Extract the data
    df <- .

    # What chromosome are we dealing with?
    chrom <- unique(df$chrom)
    # What are the traits and coefficients?
    trait <- str_c(unique(df$trait), "_", unique(df$coef))

    # Get the appropriate relationship matrix
    K_use <- K_c[[chrom]]

    # Subset the data.frame and create a model.frame
    mf <- pheno_use %>%
      select(line_name, trait) %>%
      model.frame(formula = as.formula(str_c(trait, " ~ line_name")))

    # Response vector
    y <- model.response(mf)
    # Fixed effect mean matrix
    X_mu <- model.matrix(~ 1, mf)
    # Random effect of each entry
    Z_g <- model.matrix(~ -1 + line_name, mf)
    # Fixed effect of each SNP
    X_snp <- {Z_g %*% S2TP_imputed_multi_genos_mat[,df$marker, drop = FALSE]} %>%
      remove_colinear()

    # Combine the fixed effect matrices
    X <- cbind(X_mu, X_snp)

    # Fit the model and return estimates
    fit_out <- mlmm(y = y, X = X, Z = Z_g, K = K_use)

    # Are any qvalues greater than the significance threshold?
    nonsig <- fit_out$q_value > sig_cutoff

    # While loop to backwards eliminate
    while(any(nonsig)) {

      # Which markers had the largest p_value
      which_remove <- which.max(fit_out$p_value)
      # Remove that marker from the X matrix
      X_snp <- X_snp[,-which_remove, drop = FALSE]

      # If no SNPs are present, return NA
      if (ncol(X_snp) == 0) {
        fit_out[1,] <- NA
        # Stop the loop
        break

      } else {
        # New X matrix
        X_new <- cbind(X_mu, X_snp)

        # Fit the model and return estimates
        fit_out <- mlmm(y = y, X = X_new, Z = Z_g, K = K_use)

      }

      # Are any qvalues greater than the significance threshold?
      nonsig <- fit_out$q_value > sig_cutoff

    } # End of while loop

    # Return the data.frame
    df %>%
      select(marker:pos) %>%
      right_join(., fit_out, by = "marker")
  })


# Create a K matrix using all markers
K <- A.mat(X = S2TP_imputed_multi_genos_mat, min.MAF = 0, max.missing = 1)

## Now fit a mixed model using the SNPs identified in the chromosome-specific analysis
gwas_mlmm_final_model <- gwas_mlmm_Gmodel %>%
  group_by(trait, coef) %>%
  do({
    # Extract the data
    df <- .
    
    # What are the traits and coefficients?
    trait <- str_c(unique(df$trait), "_", unique(df$coef))
    
    # Subset the data.frame and create a model.frame
    mf <- pheno_use %>%
      select(line_name, trait) %>%
      model.frame(formula = as.formula(str_c(trait, " ~ line_name")))
    
    # Response vector
    y <- model.response(mf)
    # Fixed effect mean matrix
    X_mu <- model.matrix(~ 1, mf)
    # Random effect of each entry
    Z_g <- model.matrix(~ -1 + line_name, mf)
    # Fixed effect of each SNP
    X_snp <- {Z_g %*% S2TP_imputed_multi_genos_mat[,df$marker, drop = FALSE]} %>%
      remove_colinear()
    
    # Combine the fixed effect matrices
    X <- cbind(X_mu, X_snp)
    
    # Fit the model and return estimates
    fit_out <- mlmm(y = y, X = X, Z = Z_g, K = K)
  
    # Return the data.frame
    df %>%
      select(marker:pos) %>%
      right_join(., fit_out, by = "marker")
  })


# Save this data
save_file <- file.path(result_dir, "S2MET_pheno_fw_mean_gwas_results.RData")
save("gwas_pheno_mean_fw", "gwas_mlmm_Gmodel", file = save_file)





### Mapping using resampling estimates of stability ###


## Use the FW resampling data to determine the robustness of the mapping results
# First convert to usable phenotype data.frames and output a list
resample_phenos_use <- S2MET_pheno_sample_fw %>%
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
    # Remove inifite values
    pheno_use <- core_df$pheno_use[[i]] %>% 
      mutate_at(vars(-line_name), parse_number, na = c("", NA, Inf, -Inf))
    
    gwas_out <- gwas(pheno = pheno_use, geno = genos_use, model = "G",
                     n.PC = 0, P3D = TRUE, n.core = 1)
    
    # Find the significant hits and return
    sig_gwas_list[[i]] <- gwas_out$scores %>% 
      filter(term == "main_effect") %>% 
      group_by(trait) %>% 
      mutate(q_value = qvalue(p = p_value)$qvalue) %>% 
      filter(q_value <= 0.1)
    
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
