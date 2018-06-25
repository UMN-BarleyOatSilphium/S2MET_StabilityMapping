## Variance explained by different marker subsets
## 
## Use different marker subsets to determine the variance in genotype mean and stability
## estimates that is explained by those markers
## 
## Author: Jeff Neyhart
## Last updated: June 23, 2018
## 

# # Run the source script - MSI
# repo_dir <- "/panfs/roc/groups/6/smithkp/neyha001/QTLMapping/S2MET_Mapping/"
# source(file.path(repo_dir, "source_MSI.R"))

# Run the source script - local
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

# Other packages
library(rrBLUP)
library(sommer)

# Load the FW results and marker subsets
load(file.path(result_dir, "pheno_mean_fw_results.RData"))
load(file.path(result_dir, "marker_subsets.RData"))
load(file.path(result_dir, "marker_mean_fw_results.RData"))


# Create a genotype marker matrix
M <- S2TP_imputed_multi_genos_mat[tp_geno,]
# Calculate the relationship matrix using all markers
K <- A.mat(X = M, min.MAF = 0, max.missing = 1)

# Number of model fittings
n_iter <- 100
# Empirical threshold
alpha <- 0.05

## Tidy the phenotype data for modeling
pheno_mean_fw_tidy <- pheno_mean_fw %>% 
  distinct(trait, line_name, g, b, delta) %>% 
  mutate(log_delta = log(delta)) %>% 
  dplyr::select(-delta) %>%
  gather(coef, value,  g:log_delta)



##
## First use all markers to calculate variance components
##

pheno_fw_var_comp_all_markers <- pheno_mean_fw_tidy %>% 
  group_by(trait, coef) %>%
  do({
    df <- .
    Z <- model.matrix(~ -1 + line_name, df)
    fit <- mixed.solve(y = df$value, Z = Z, K = K)
    
    data.frame(snps = fit$Vu, res = fit$Ve)
    
  }) %>% ungroup()




##
## Second use the 95% threhold to determine plastic and stable markers
## Run variance component analysis using random subsets of the stable markers
##

# Find the most plastic markers and stable markers
marker_mean_fw_tidy <- marker_mean_fw %>% 
  distinct(marker, chrom, pos, trait, c, delta) %>%
  dplyr::select(-delta)


# For each trait, calculate empirical thresholds for significance
marker_fw_sig <- marker_mean_fw_tidy %>%
  group_by(trait) %>% 
  mutate(lower_perc = quantile(c, alpha / 2), 
         upper_perc = quantile(c, 1 - (alpha / 2))) %>%
  ungroup() %>%
  mutate(significance = case_when(c >= upper_perc ~ "plastic",
                                  c <= lower_perc ~ "plastic",
                                  TRUE ~ "stable"),
         marker_type = if_else(str_detect(marker, "^S"), "GBS", "BOPA"))




# Create a list
var_comp_list <- list()

## Iterate over the traits
for (tr in unique(marker_fw_sig$trait)) {

  stab_mar <- marker_fw_sig %>%
    subset(significance == "stable" & trait == tr, marker, drop = TRUE)
  
  plas_mar <- marker_fw_sig %>%
    subset(significance == "plastic" & trait == tr, marker, drop = TRUE)
  
  # What is the minimum number of markers to sample?
  min_marker <- min(c(length(stab_mar), length(plas_mar)))

  stab_mar_samples <- replicate(n = n_iter, sample(x = stab_mar, size = min_marker), simplify = FALSE)

  ## Subset marker matrices and create the static relationship matrices
  K_plas <- M[,plas_mar,drop = FALSE] %>% A.mat(X = ., min.MAF = 0, max.missing = 1)

  # List of non-sig K matrices
  K_stab_sample <- stab_mar_samples %>%
    map(~M[,.,drop = FALSE]) %>%
    map(~A.mat(X = ., min.MAF = 0, max.missing = 1))
  
  
  ## Subset the data to model
  pheno_mean_fw_tidy_tomodel <- pheno_mean_fw_tidy %>%
    filter(trait == tr) %>%
    split(.$coef)
  
  var_comp_out <- vector("list", n_iter)
    
  
  ## Model
  var_comp_out <- K_stab_sample %>%
    lapply(function(K_stab) {
      
      # For each trait-coef, create mf, model.matrices, and fit
      pheno_mean_fw_tidy_tomodel %>% 
        map(~{
          mf <- model.frame(value ~ line_name, .)
          y <- model.response(mf)
          Z <- model.matrix(~ -1 + line_name, mf) %>%
            `colnames<-`(., row.names(M))
          X <- model.matrix(~ 1, mf)
          
          fit <- mmer(Y = y, X = X, Z = list(stab = list(Z = Z, K = K_stab), plas = list(Z = Z, K = K_plas)), silent = TRUE)
          
          fit$var.comp %>% map(as.numeric) %>% as.data.frame() })

    })
    
  # Combine into a data.frame, then add to a list
  var_comp_list[[tr]] <- var_comp_out %>% 
    transpose() %>%
    map(~bind_rows(.) %>% mutate(iter = seq(nrow(.)))) %>%
    list(., names(.)) %>%
    pmap_df(~mutate(.x, coef = .y)) %>%
    mutate(trait = tr) %>%
    dplyr::select(trait, coef, iter, names(.))
  
}

# Tidy up
pheno_fw_var_comp_samples <- var_comp_list %>% 
  bind_rows()


# Save the results
save("pheno_fw_var_comp_all_markers", "pheno_fw_var_comp_samples", 
     file = file.path(result_dir, str_c("pheno_mar_fw_varcomp.RData")))
