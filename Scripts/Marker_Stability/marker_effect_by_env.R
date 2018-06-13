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

n_core <- 7

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


snps_by_chrom <- map(snps_by_chrom, head, 10)


## Split SNPs by chromosome
K_chr <- snps_by_chrom %>%
  map(~setdiff(snp_info$marker, .)) %>% 
  map(~A.mat(X = M[,.], min.MAF = 0, max.missing = 1))

# Create an empty list for storing marker coefficients
marker_by_env_effects <- list()

## Fit the mixed model to estimate variance components
# Iterate over traits
for (tr in unique(S2_MET_BLUEs_use$trait)) {
  # Print the trait
  print(tr)
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
  
  # Create a large list to iterate
  parallel_list <- transpose(list(snps_by_chrom, K_chr, population_parameters))
  
  # Parallelize
  marker_score <- mclapply(X = parallel_list, FUN = function(chrom_list) {
    
      # Population structure matrix
      # 1 PC
      Q <- eigen(chrom_list[[2]])$vectors[,1]
      X1 <- cbind(X, pop_str = (Z %*% Q))
      Hinv <- chrom_list[[3]]
      
      apply(X = M[,chrom_list[[1]]], MARGIN = 2, FUN = function(snp) {
        
        # Create the SNP X E matrix
        X_snp <- X_env * c(Z %*% snp)

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
      
    }, mc.cores = n_core)
  
  # # Rotate and convert to df
  # marker_by_env_effects[[tr]] <- marker_score %>% 
  #   map(~as.data.frame(t(.)) %>% rownames_to_column("marker")) %>% 
  #   bind_rows() %>%
  #   rename_at(vars(-marker), ~str_extract(., "[A-Z]{3}[0-9]{2}")) %>%
  #   gather(environment, effect, -marker) %>%
  #   mutate(trait = tr)
  
  # Return the matrix
  marker_by_env_effects[[tr]] <- marker_score
  
}

## Save the population parameters
save_file <- file.path(result_dir, "marker_by_env_effects.RData")
save("marker_by_env_effects", file = save_file)



# ## Calculate marker effects separately for each environment
# marker_by_env <- S2_MET_BLUEs_use %>%
#   group_by(trait, environment) %>%
#   do({
#     df <- .
#     
#     # Model frame
#     mf <- model.frame(value ~ line_name + std_error, df)
#     
#     # Response
#     y <- model.response(mf)
#     # Grand mean design
#     X <- model.matrix(~ 1, mf)
#     # Line_name design
#     Zg <- model.matrix(~ -1 + line_name, mf)
#     
#     # Subset the marker matrix
#     Z <- Zg %*% M
#     K <- diag(ncol(Z))
#     
#     # Pull out the weights into an R matrix
#     # R <- solve(diag(mf$std_error^2))
#     
#     # Fit
#     # fit <- sommer::mmer(Y = y, X = X, Z = list(marker = list(Z = Z, K = K)), R = list(res = R))
#     fit <- mixed.solve(y = y, Z = Z, method = "REML")
#     
#     
#     # # Grab the marker effects and return
#     # marker_effect <- fit$u.hat$marker %>% 
#     #   as.data.frame() %>% 
#     #   rownames_to_column("marker") %>% 
#     #   dplyr::select(marker, effect = T1)
#     
#     # Grab the marker effects and return
#     marker_effect <- fit$u %>% 
#       as.data.frame() %>% 
#       rownames_to_column("marker") %>%
#       rename(effect = ".")
#     
#     # Return a data_frame
#     data_frame(marker_effect = list(marker_effect), fit = list(fit)) })


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


