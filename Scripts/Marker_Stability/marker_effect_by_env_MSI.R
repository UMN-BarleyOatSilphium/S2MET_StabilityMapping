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


## Fit the mixed model to estimate variance components
population_parameters <- list()

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
  population_parameters[[tr]] <- list()
  
  for (i in seq_along(K_chr)) {
    Q <- eigen(K_chr[[i]])$vectors[,1]
    X1 <- cbind(X, Z %*% Q)
    population_parameters[[tr]][[as.character(i)]] <- 
      mixed.solve(y = y, Z = Z, K = K_chr[[i]], X = X1, return.Hinv = TRUE)$Hinv
    
  }
  
}


## Save the population parameters
save_file <- file.path(result_dir, "marker_by_env_population_parameters.RData")
save("population_parameters", file = save_file)

