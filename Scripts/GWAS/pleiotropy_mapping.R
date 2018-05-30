## S2MET GWAS for pleiotropy between genotype mean and stability
## 
## Perform pleiotropic association analysis of the FW stability coefficients and 
## genotype means
## 
## Author: Jeff Neyhart
## Last updated: May 30, 2018
## 


# Run the source script - MSI
repo_dir <- "/panfs/roc/groups/6/smithkp/neyha001/QTLMapping/S2MET_Mapping/"
source(file.path(repo_dir, "source_MSI.R"))

# # Run the source script - local
# repo_dir <- getwd()
# source(file.path(repo_dir, "source.R"))


# Extract the SNP information for later use
snp_info <- S2TP_imputed_multi_genos_hmp %>%
  select(marker = rs, chrom, pos, cM_pos)


# Load the FW results
load(file.path(result_dir, "S2MET_pheno_mean_fw_results.RData"))
# Load the FW sampling results
load(file.path(result_dir, "S2MET_pheno_fw_resampling.RData"))

## Number of cores
n_core <- ifelse(Sys.info()["sysname"] == "Windows", 1, detectCores())


# Format the genotype data for use
geno_use <- S2TP_imputed_multi_genos_hmp %>%
  select(-alleles, -cM_pos) %>%
  rename(marker = rs) %>%
  as.data.frame()


# Format for modeling
pheno_to_model <- S2_MET_pheno_mean_fw %>% 
  distinct(line_name, trait, g, b, delta) %>% 
  mutate(log_delta = log(delta)) %>%
  select(-delta) %>%
  split(.$trait) %>%
  map(select, -trait) %>%
  map(as.data.frame)

# K matrix
M <- S2TP_imputed_multi_genos_mat
K <- A.mat(X = M, min.MAF = 0, max.missing = 1)
# Q vector
Q <- eigen(x = K)$vector[,1]

# K matrix by chromosome
markers_per_chrom <- geno_use %>% split(.$chrom) %>% map("marker")
K_per_chrom <- markers_per_chrom %>% map(~setdiff(geno_use$marker, .)) %>% 
  map(~A.mat(X = S2TP_imputed_multi_genos_mat[,.], min.MAF = 0, max.missing = 1))
Q_per_chrom <- map(K_per_chrom, ~eigen(x = .)$vector[,1])



## Models to fit
# models <- c("K", "QK", "G", "QG")
models <- "QK"


## Test for pleiotropy
##
## # Fit a multi-variate mixed model for the genotype mean and stability,
## # then test each marker using an intersection-union test
##
## # H0: the marker is not associated with the mean or it is not associated with stability
## # HA: the marker is associated with both mean and stability
##

# Re-organize the data for modeling
pheno_to_model_plei <- pheno_to_model %>% 
  list(., names(.)) %>% 
  pmap_df(~mutate(., trait = .y)) %>% 
  gather(stab_coef, value, b, log_delta)

# Split up the data for iterations
pheno_to_model_plei_split <- pheno_to_model_plei %>%
  split(list(.$trait, .$stab_coef))

# Create an empty list to store results
gwas_pheno_mean_fw_plei_out <- vector("list", length = length(pheno_to_model_plei_split)) %>%
  set_names(., names(pheno_to_model_plei_split))


# Split the chromosomes by the number of cores
snps_by_core <- snp_info %>%
  assign_cores(n_core) %>%
  split(.$core)



# Iterate over trait-stability combinations
for (i in seq_along(pheno_to_model_plei_split)) {
  
  # Extract the data
  df <- pheno_to_model_plei_split[[i]]
  
  # For each trait and parameter, conduct an association test for each marker
  # Create a matrix of Y
  Y <- df %>%
    select(g, value) %>%
    as.matrix() %>%
    t() # Transpose for the function
  
  X_mu <- model.matrix(~ 1, df)
  Z <- t(model.matrix(~ -1 + line_name, df))
  
  # Iterate over the markers split up by core
  core_out <- mclapply(X = snps_by_core, FUN = function(marker_core) {
    
    # Empty list for output
    marker_core_out <- vector("list", length = nrow(marker_core))
    # Add models
    marker_core1 <- crossing(marker_core, model = models)
    
    # Iterate
    for (j in seq(nrow(marker_core1))) {
      
      # Add the Q vector?
      if (grepl(pattern = "Q", x = marker_core1$model[j])) {
        q <- Q
      } else {
        q <- NULL
      }
      
      # Combine the snp with the mean
      X <- cbind(X_mu, M[,marker_core1$marker[j],drop = FALSE], Q = q)
      
      # Fit
      fit <- emmremlMultivariate(Y = Y, X = t(X), Z = Z, K = K, varBhat = T)
      marker_core_out[[j]] <- data.frame(marker = marker_core1$marker[j],
                                         term = row.names(fit$Bhat), beta = fit$Bhat[,2], se = sqrt(fit$varBhat[3:4]),
                                         row.names = NULL, stringsAsFactors = FALSE)
      
    }
    
    # Return the output
    return(marker_core_out)
    
  }, mc.cores = n_core)
  
  # Add the results to the list
  gwas_pheno_mean_fw_plei_out[[i]] <- map_df(core_out, bind_rows)
  
}

## Save the data
save_file <- file.path(result_dir, "pheno_fw_gwas_pleiotropy_results.RData")
save("gwas_pheno_mean_fw_plei_out", file = save_file)





