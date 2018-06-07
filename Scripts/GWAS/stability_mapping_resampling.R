## S2MET GWAS of Stability Coefficients - resampling
## 
## Perform association analysis of the FW stability coefficients for regular phenotype
## FW, and for FW on environmental covariables
## 
## Author: Jeff Neyhart
## Last updated: June 01, 2018
## 


# Run the source script - MSI
repo_dir <- "/panfs/roc/groups/6/smithkp/neyha001/QTLMapping/S2MET_Mapping/"
source(file.path(repo_dir, "source_MSI.R"))

# # Run the source script - local
# repo_dir <- getwd()
# source(file.path(repo_dir, "source.R"))


# Load the FW results
load(file.path(result_dir, "pheno_mean_fw_results.RData"))
# Load the FW sampling results
load(file.path(result_dir, "pheno_fw_resampling.RData"))


## Number of cores
n_core <- ifelse(Sys.info()["sysname"] == "Windows", 1, detectCores())

alpha <- 0.05

# Format the genotype data for use
geno_use <- S2TP_imputed_multi_genos_hmp %>%
  select(-alleles, -cM_pos) %>%
  rename(marker = rs) %>%
  as.data.frame()


## Use the FW resampling data to determine the robustness of the mapping results
# Tidy up for splitting by cores
resample_phenos_use <- pheno_sample_mean_fw %>%
  # Convert delta to log_delta
  mutate(log_delta = log(delta)) %>% 
  # Tidy
  select(-delta) %>% 
  group_by(trait, p, iter) %>% 
  nest()

# Assign cores
resample_phenos_use_list <- resample_phenos_use %>%
  assign_cores(n_core = n_core) %>%
  split(.$core)



## Parallelize the association mapping runs using the resampling
resample_gwas_sig_out <- mclapply(X = resample_phenos_use_list, function(core_df) {
  
  # Output list
  results_out <- vector("list", nrow(core_df))
  
  # Iterate over rows
  for (i in seq(nrow(core_df))) {
    
    pheno <- as.data.frame(core_df$data[[i]])
    
    # Run the GWAS
    gwas_out <- GWAS(pheno = pheno, geno = geno_use, n.PC = 1, plot = FALSE)
    
    # Grab the scores, adjust the pvalues, then filter for those <= alpha
    gwas_out_tidy <- gwas_out %>% 
      gather(coef, neg_log_p, b, log_delta) %>% 
      mutate(pvalue = 10^-neg_log_p)
    
    # Filter and adjust pvalues
    gwas_sig_scores <- gwas_out_tidy %>%
      filter(pvalue < 1) %>% # Remove the p == 1 (missing)
      mutate(padj = p.adjust(pvalue, method = "fdr")) %>%
      filter(padj <= alpha)
    
    # Return a list
    results_out[[i]] <- list(gwas_sig = gwas_sig_scores, n_NA = group_by(gwas_out_tidy, coef) %>% 
                               summarize(n_NA = sum(pvalue == 1)))
    
  }
  
  # Add the results to the core_df and return
  core_df %>% 
    mutate(results = results_out) %>%
    select(-data, -core)
  
}, mc.cores = n_core)



# Save the results
save_file <- file.path(result_dir, str_c("pheno_fw_gwas_resample_results.RData"))
save("resample_gwas_sig_out", file = save_file)

