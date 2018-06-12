## S2MET Mapping
## Calculate marker effects in each environment
## 
## This script will use MSI to calculate the Hinv matrix for the multi-environment
## data
## 
## Author: Jeff Neyhart
## Last updated: June 11, 2018
## 


# Run the source script - MSI
repo_dir <- "/panfs/roc/groups/6/smithkp/neyha001/QTLMapping/S2MET_Mapping/"
source(file.path(repo_dir, "source_MSI.R"))

# # Run the source script - local
# repo_dir <- getwd()
# source(file.path(repo_dir, "source.R"))


# Rename the marker matrix
M <- S2TP_imputed_multi_genos_mat

# Get the marker information
snp_info <- S2TP_imputed_multi_genos_hmp %>%
  select(marker = rs, chrom, pos, cM_pos) %>%
  # Correct the cM position for BOPA snps
  mutate(cM_pos = if_else(str_detect(marker, "^S"), cM_pos, cM_pos / 1000))

# Subset the entry list for the tp
tp_entry_list <- entry_list %>%
  filter(Line %in% tp_geno) %>%
  select(line_name = Line, program = Program)

# Significance threshold for GWAS
alpha <- 0.05

# What model will we use for mapping?
model_use <- "QG"


## Use the GWAS Q+G model to estimate the effect of each marker in each environment
## Use P3D

## Split SNPs by chromosomes
snps_by_chrom <- snp_info %>% 
  split(.$chrom) %>%
  map("marker")


## Split SNPs by chromosome
K_chr <- snps_by_chrom %>%
  map(~setdiff(snp_info$marker, .)) %>% 
  map(~A.mat(X = M[,.], min.MAF = 0, max.missing = 1))

# Create an empty list for storing marker coefficients
marker_by_env_effects <- list()

## Fit the mixed model to estimate variance components
# Iterate over traits
for (tr in unique(S2_MET_BLUEs_use$trait)) {
  
  df <- subset(S2_MET_BLUEs_use, trait == tr)
  
  ## May need to edit this with the correct weights
  mf <- model.frame(value ~ line_name + environment, df)
  y <- model.response(mf)
  
  # Fixed effects of environment
  X <- model.matrix(~ environment, mf)
  Z <- model.matrix(~ -1 + line_name, mf)
  
  # Iterate over the K matrices and return the H_inv matrix
  population_parameters <- list()
  
  for (i in seq_along(K_chr)) {
    Q <- eigen(K_chr[[i]])$vectors[,1]
    X1 <- cbind(X, Z %*% Q)
    population_parameters[[as.character(i)]] <- 
      mixed.solve(y = y, Z = Z, K = K_chr[[i]], X = X1, return.Hinv = TRUE)$Hinv
    
  }
  
  # Create the snp:environment model matrix
  X_env <- model.matrix(~ -1 + environment, mf)
  
  ## Use the GWAS QG model to estimate the effect of each marker in each environment
  ## This will be the environment-specific marker effect + the mean
  
  snps_by_chrom1 <- snps_by_chrom %>% map(head)
  
  # Iterate over the snps per chromosome
  marker_score <- list(snps_by_chrom1, K_chr, population_parameters) %>%
    pmap(~{
      
      # Population structure matrix
      # 1 PC
      Q <- eigen(..2)$vectors[,1]
      X1 <- cbind(X, pop_str = (Z %*% Q))
      Hinv <- ..3
      
      apply(X = M[,..1], MARGIN = 2, FUN = function(snp) {
        
        # Create the SNP X E matrix
        X_snp <- X_env * c(Z %*% snp)
        colnames(X_snp) <- paste0(colnames(X_snp), "_snp")
        
        # Create a new X matrix
        X2 <- cbind(X1, X_snp)
        
        # Index of snp betas
        j <- seq(ncol(X1) + 1, ncol(X2))
        
        W <- crossprod(X2, Hinv %*% X2)
        Winv <- try(solve(W), silent = TRUE)
        
        
        if (class(Winv) != "try-error") {
          beta <- Winv %*% crossprod(X2, Hinv %*% y)
          beta <- beta[j,]
          
        } else {
          beta <- NA
          
        }
        
        return(beta) })
      
    })
  
  # Rotate and convert to df
  marker_by_env_effects[[tr]] <- marker_score %>% 
    map(~as.data.frame(t(.)) %>% rownames_to_column("marker")) %>% 
    bind_rows() %>%
    rename_at(vars(-marker), ~str_extract(., "[A-Z]{3}[0-9]{2}")) %>%
    gather(environment, effect, -marker)
  
}

## Save the population parameters
save_file <- file.path(result_dir, "marker_by_env_effects.RData")
save("marker_by_env_effects", file = save_file)

