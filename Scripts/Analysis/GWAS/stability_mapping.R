## S2MET GWAS of Stability Coefficients
## 
## Perform association analysis of the FW stability coefficients for regular phenotype
## FW, and for FW on environmental covariables
## 
## 


# List of packages to load
packages <- c("dplyr", "purrr", "tibble", "tidyr", "readr", "stringr", "readxl", "modelr", 
              "parallel", "purrrlyr", "rrBLUP", "EMMREML", "qvalue")


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

geno_dir <- bopa_geno_dir <-  "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/Data/Genos"
pheno_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/Data/Phenos"


# Other directories
script_dir <- file.path(proj_dir, "Scripts/")
analysis_dir <- file.path(script_dir, "Analysis")
fig_dir <- file.path(script_dir, "Figures/")
map_dir <- file.path(analysis_dir, "GWAS")
result_dir <- file.path(proj_dir, "Results")


# # Load the phenotypic data
# load(file.path(pheno_dir, "S2_MET_BLUEs.RData"))
# Load the genotypic data
load(file.path(bopa_geno_dir, "S2TP_multi_genos.RData"))

# Extract the SNP information for later use
snp_info <- S2TP_imputed_multi_genos_hmp %>%
  select(marker = rs, chrom, pos, cM_pos)


# Load the FW results
load(file.path(result_dir, "S2MET_pheno_mean_fw_results.RData"))
# Load the FW sampling results
load(file.path(result_dir, "S2MET_pheno_fw_resampling.RData"))

# Number of cores
n_core <- detectCores()



# Format for modeling
pheno_to_model <- S2_MET_pheno_mean_fw %>% 
  distinct(line_name, trait, g, b, delta) %>% 
  mutate(log_delta = log(delta)) %>%
  select(-delta) %>%
  gather(coef, value, g, b, log_delta) 


## Association analysis
# Fit model chromosome-wise (G model)
snps_by_chrom <- snp_info %>%
  split(.$chrom) %>%
  map("marker")

# Marker matrices per chromosome
M_chr <- snps_by_chrom %>% 
  map(~S2TP_imputed_multi_genos_mat[,., drop = FALSE])

# Create relationship matrices for each chromosome
K_chr <- snps_by_chrom %>% 
  map(~S2TP_imputed_multi_genos_mat[,setdiff(colnames(S2TP_imputed_multi_genos_mat), .), drop = FALSE]) %>%
  map(., ~A.mat(., min.MAF = 0, max.missing = 1))


# # For each trait and parameter, conduct an association test for each marker
# gwas_out <- pheno_to_model %>%
#   group_by(trait, coef) %>%
#   do({
#     df <- .
# 
#     # EMMREML
#     y <- df$value
#     X_mu <- model.matrix(~ 1, df)
#     Z <- model.matrix(~ -1 + line_name, df)
#     
#     # Iterate over the marker list and the K matrix list
#     fit_out <- pmap(list(M_chr, K_chr), function(.x, .y) {
#       # Apply over markers
#       apply(X = .x, MARGIN = 2, FUN = function(snp) {
#         
#         # Combine the snp with the mean
#         X <- cbind(X_mu, snp)
#         
#         # Fit and return
#         fit <- emmreml(y = y, X = X, Z = Z, K = .y, varbetahat = T)
#         fit[c("Vu", "Ve", "betahat", "varbetahat", "loglik")]
#         
#       }) })
#     
#     # Gather the marker betas and standard errors
#     mar_effect_estimates <- fit_out %>%
#       map(~map(., ~data.frame(beta = .$betahat[2], se = sqrt(.$varbetahat[2]), row.names = NULL))) %>%
#       map_df(~pmap_df(list(., names(.)), ~mutate(.x, marker = .y)))
#     
#     # Return
#     mar_effect_estimates %>%
#       select(marker, names(.)) })
# 
# # Calculate test statistics and pvalues
# gwas_pheno_mean_fw <- gwas_out %>%
#   ungroup() %>%
#   mutate(statistic = (beta^2) / (se^2), 
#          pvalue = pchisq(q = statistic, df = 1, lower.tail = FALSE)) %>%
#   # Add the SNP metadata
#   full_join(., select(S2TP_imputed_multi_genos_hmp, marker = rs, alleles:cM_pos)) %>%
#   select(marker, alleles:cM_pos, trait, coef, beta:pvalue)
#     
# 
# 
# ## Fit a multi-locus mixed model to identify significant QTL
# 
# # First correct the p_values
# gwas_pheno_mean_fw_adj <- gwas_pheno_mean_fw %>%
#   group_by(trait, coef) %>%
#   mutate(p_adj = p.adjust(pvalue, method = "fdr"),
#          q_value = qvalue(p = pvalue)$qvalue,
#          local_fdr = qvalue(p = pvalue)$lfdr,
#          neg_log10_fdr05 = -log10(0.05),
#          neg_log10_fdr10 = -log10(0.10)) %>%
#   mutate_at(vars(pvalue, p_adj, q_value), funs(neg_log10 = -log10(.))) %>%
#   ungroup()
# 
# # Cutoff for significance
# sig_cutoff <- 0.05
# 
# # Select the G model and identify the signficant associations
# gwas_adj_sig <- gwas_pheno_mean_fw_adj %>%
#   # filter(model == "G", q_value <= sig_cutoff)
#   filter(q_value <= sig_cutoff)
# 
# 
# 
# ## Fit MLMM on a chromosome-wise basis
# # Define a function that fits the linear mm and returns estimates, std.errors, and pvalue of
# # markers
# mlmm <- function(y, X, Z, K) {
#   # Fit the model
#   fit <- mixed.solve(y = y, Z = Z, K = K, X = X, SE = TRUE)
#   alpha <- fit$beta[-1]
#   se <- fit$beta.SE[-1]
#   # Perform hypothesis test
#   chisq <- alpha^2 / se^2
#   p_value <- pchisq(q = chisq, df = 1, lower.tail = FALSE)
#   q_value <- p.adjust(p = p_value, method = "BH")
# 
#   # Create an empty vector
#   R_sqr_ind <- numeric(ncol(X) - 1)
#   
#   for (i in seq_along(R_sqr_ind)) {
#     X_use <- X[,c(1, i + 1), drop = FALSE]
#     y_hat <- (X_use %*% fit$beta[c(1, i + 1)])
#     R_sqr_ind[i] <- cor(y, y_hat)^2
#   }
#   
#   y_hat_snp <- (X %*% fit$beta)
#   
#   # Calculate R^2 and adjusted R^2 for the full model
#   n <- length(y)
#   p <- ncol(X) - 2
# 
#   R_sqr <- cor(y, y_hat_snp)^2
#   R_sqr_adj <- R_sqr - ((1 - R_sqr) * (p / (n - p - 1)))
# 
#   # Return a list and a data.frame
#   data.frame(marker = names(alpha), alpha, se, p_value, q_value, R_sqr_snp = R_sqr_ind,
#              R_sqr = R_sqr, R_sqr_adj = R_sqr_adj, stringsAsFactors = FALSE)
# 
# }
# 
# # Function for removing colinear SNPs
# remove_colinear <- function(x) {
#   # LD
#   x_LD <- cor(x)^2
#   # Find the index of the SNPs to be dropped
#   snp_index <- apply(X = x_LD, MARGIN = 1, FUN = function(LD) max(which(LD > 0.99)))
#   to_drop <- which(duplicated(snp_index))
# 
#   ifelse(length(to_drop) > 0, return(x[,-to_drop, drop = FALSE]), return(x))
# 
# }
# 
# # Iterate over trait, coef, and chromosome
# gwas_mlmm_Gmodel <- gwas_adj_sig %>%
#   group_by(trait, coef, chrom) %>%
#   do({
#     # Extract the data
#     df <- .
# 
#     # What chromosome are we dealing with?
#     chr <- unique(df$chrom)
#     # What are the traits and coefficients?
#     tr <- unique(df$trait)
#     cf <- unique(df$coef)
#     
#     # Get the appropriate relationship matrix
#     K_use <- K_chr[[chrom]]
# 
#     # Subset the data.frame and create a model.frame
#     mf <- pheno_to_model %>%
#       filter(trait == tr, coef == cf) %>%
#       model.frame(formula = value ~ line_name)
# 
#     # Response vector
#     y <- model.response(mf)
#     # Fixed effect mean matrix
#     X_mu <- model.matrix(~ 1, mf)
#     # Random effect of each entry
#     Z_g <- model.matrix(~ -1 + line_name, mf)
#     # Fixed effect of each SNP
#     X_snp <- {Z_g %*% M[,df$marker, drop = FALSE]} %>%
#       remove_colinear()
# 
#     # Combine the fixed effect matrices
#     X <- cbind(X_mu, X_snp)
# 
#     # Fit the model and return estimates
#     fit_out <- mlmm(y = y, X = X, Z = Z_g, K = K_use)
# 
#     # Are any qvalues greater than the significance threshold?
#     nonsig <- fit_out$q_value > sig_cutoff
# 
#     # While loop to backwards eliminate
#     while(any(nonsig)) {
# 
#       # Which markers had the largest p_value
#       which_remove <- which.max(fit_out$p_value)
#       # Remove that marker from the X matrix
#       X_snp <- X_snp[,-which_remove, drop = FALSE]
# 
#       # If no SNPs are present, return NA
#       if (ncol(X_snp) == 0) {
#         fit_out[1,] <- NA
#         # Stop the loop
#         break
# 
#       } else {
#         # New X matrix
#         X_new <- cbind(X_mu, X_snp)
# 
#         # Fit the model and return estimates
#         fit_out <- mlmm(y = y, X = X_new, Z = Z_g, K = K_use)
# 
#       }
# 
#       # Are any qvalues greater than the significance threshold?
#       nonsig <- fit_out$q_value > sig_cutoff
# 
#     } # End of while loop
# 
#     # Return the data.frame
#     df %>%
#       select(marker:pos) %>%
#       right_join(., fit_out, by = "marker")
#   })
# 
# 
# # Create a K matrix using all markers
# K <- A.mat(X = S2TP_imputed_multi_genos_mat, min.MAF = 0, max.missing = 1)
# 
# ## Now fit a mixed model using the SNPs identified in the chromosome-specific analysis
# gwas_mlmm_final_model <- gwas_mlmm_Gmodel %>%
#   group_by(trait, coef) %>%
#   do({
#     # Extract the data
#     df <- .
#     
#     # What are the traits and coefficients?
#     tr <- unique(df$trait)
#     cf <- unique(df$coef)
#     
#     # Subset the data.frame and create a model.frame
#     mf <- pheno_to_model %>%
#       filter(trait == tr, coef == cf) %>%
#       model.frame(formula = value ~ line_name)
#     
#     # Response vector
#     y <- model.response(mf)
#     # Fixed effect mean matrix
#     X_mu <- model.matrix(~ 1, mf)
#     # Random effect of each entry
#     Z_g <- model.matrix(~ -1 + line_name, mf)
#     # Fixed effect of each SNP
#     X_snp <- {Z_g %*% S2TP_imputed_multi_genos_mat[,df$marker, drop = FALSE]} %>%
#       remove_colinear()
#     
#     # Combine the fixed effect matrices
#     X <- cbind(X_mu, X_snp)
#     
#     # Fit the model and return estimates
#     fit_out <- mlmm(y = y, X = X, Z = Z_g, K = K)
#   
#     # Return the data.frame
#     df %>%
#       select(marker:pos) %>%
#       right_join(., fit_out, by = "marker")
#   })
# 
# 
# # Save this data
# save_file <- file.path(result_dir, "S2MET_pheno_fw_mean_gwas_results.RData")
# save("gwas_pheno_mean_fw", "gwas_mlmm_Gmodel", "gwas_mlmm_final_model", file = save_file)
# 
# 
# 
# 
# 
### Mapping using resampling estimates of stability ###


## Use the FW resampling data to determine the robustness of the mapping results
# Tidy up for splitting by cores
resample_phenos_use <- S2MET_pheno_sample_fw %>%
  # Convert delta to log_delta
  mutate(log_delta = log(delta)) %>% 
  # Tidy
  select(-delta) %>% 
  gather(term, estimate, b, log_delta) %>% 
  group_by(trait, p, iter, term) %>% 
  nest()

# Assign cores
resample_phenos_use_list <- resample_phenos_use %>%
  ungroup() %>%
  mutate(core = sort(rep(seq(n_core), length.out = nrow(.)))) %>%
  split(.$core)

# Write a function to conduct the association analysis on each marker
marker_assoc <- function(data, M_chr, K_chr) {
 
  # Create model matrices
  y <- data$estimate
  X_mu <- model.matrix(~ 1, data)
  Z <- model.matrix(~ -1 + line_name, data)
  
  # Iterate over the marker list and the K matrix list
  fit_out <- pmap(list(M_chr, K_chr), function(.x, .y) {
    # Apply over markers
    apply(X = .x, MARGIN = 2, FUN = function(snp) {
      
      # Combine the snp with the mean
      X <- cbind(X_mu, snp)
      
      # Fit and return
      fit <- emmreml(y = y, X = X, Z = Z, K = .y, varbetahat = T)
      fit[c("Vu", "Ve", "betahat", "varbetahat", "loglik")]
      
      data.frame(beta = fit$betahat[2,1], se = sqrt(fit$varbetahat[2]),
                 row.names = NULL, stringsAsFactors = FALSE)
      
    }) })
  
  # Calculate p-values
  fit_out %>% 
    map_df(~list(., names(.)) %>% pmap_df(~mutate(.x, marker = .y))) %>%
    mutate(statistic = (beta^2) / (se^2),
           pvalue = pchisq(q = statistic, df = 1, lower.tail = FALSE), 
           padj = p.adjust(p = pvalue, method = "fdr")) %>% 
    filter(padj <= 0.1)
  
}

# Iterate over the list of data.frames and perform the association
# Then filter for the significant markers

# Parallelize
resample_gwas_sig_out <- mclapply(X = resample_phenos_use_list, FUN = function(pheno_model) {
  
  # Map over the data
  data_list <- pheno_model$data
  
  gwas_out_list <- map(data_list, ~marker_assoc(data = ., M_chr = M_chr, K_chr = K_chr))
  
  # Add snp information
  gwas_out_list_ann <- map(gwas_out_list, ~right_join(snp_info, ., by = "marker"))
  
  # Add the output to the pheno_model df and return
  pheno_model %>%
    select(trait:term) %>%
    mutate(gwas_out = gwas_out_list_ann)
  
})
  


# Save the results
save_file <- file.path(result_dir, "S2MET_pheno_fw_gwas_resample_results.RData")
save("resample_gwas_sig_out", file = save_file)




# 
# ## Test for pleiotropy
# ## 
# ## # Fit a multi-variate mixed model for the genotype mean and stability,
# ## # then test each marker using an intersection-union test
# ## 
# ## # H0: the marker is not associated with the mean or it is not associated with stability
# ## # HA: the marker is associated with both mean and stability
# ## 
# 
# # Re-organize the data for modeling
# pheno_to_model_plei <- pheno_to_model %>%
#   spread(coef, value) %>% 
#   gather(stab_coef, value, b, log_delta)
# 
# # Split up the data for iterations
# pheno_to_model_plei_split <- pheno_to_model_plei %>% 
#   split(list(.$trait, .$stab_coef))
# 
# # Create an empty list to store results
# gwas_pheno_mean_fw_plei_out <- vector("list", length = length(pheno_to_model_plei_split)) %>%
#   set_names(., names(pheno_to_model_plei_split))
# 
# 
# # Split the chromosomes by the number of cores
# snps_by_core <- snp_info %>% 
#   mutate(core = sort(rep(seq(n_core), length.out = nrow(.)))) %>%
#   split(.$core)
# 
# 
# # Iterate over trait-stability combinations
# for (i in seq_along(pheno_to_model_plei_split)) {
#   
#   # Extract the data
#   df <- pheno_to_model_plei_split[[i]]
#   
#   # For each trait and parameter, conduct an association test for each marker
#   # Create a matrix of Y
#   Y <- df %>%
#     select(g, value) %>%
#     as.matrix() %>%
#     t() # Transpose for the function
#     
#   X_mu <- model.matrix(~ 1, df)
#   Z <- model.matrix(~ -1 + line_name, df) %>%
#     t()
#   
#   # Iterate over the markers split up by core
#   core_out <- mclapply(X = snps_by_core, FUN = function(marker_core) {
#     
#     # Empty list for output
#     marker_core_out <- vector("list", length = nrow(marker_core))
#     
#     # Iterate
#     for (j in seq(nrow(marker_core))) {
#       
#       # Combine the snp with the mean
#       X <- cbind(X_mu, S2TP_imputed_multi_genos_mat[,marker_core$marker[j]])
#       
#       # Fit
#       fit <- emmremlMultivariate(Y = Y, X = t(X), Z = Z, K = K_chr[[marker_core$chrom[j]]], varBhat = T)
#       marker_core_out[[j]] <- data.frame(marker = marker_core$marker[j],
#         term = row.names(fit$Bhat), beta = fit$Bhat[,2], se = sqrt(fit$varBhat[3:4]),
#         row.names = NULL, stringsAsFactors = FALSE)
#       
#     }
#     
#     # Return the output
#     return(marker_core_out)
#     
#   }, mc.cores = n_core)
#   
#   # Add the results to the list
#   gwas_pheno_mean_fw_plei_out[[i]] <- map_df(core_out, bind_rows)
#   
# }
#   
# ## Save the data
# save_file <- file.path(result_dir, "S2MET_pheno_fw_gwas_pleiotropy_results.RData")
# save("gwas_pheno_mean_fw_plei_out", file = save_file)
