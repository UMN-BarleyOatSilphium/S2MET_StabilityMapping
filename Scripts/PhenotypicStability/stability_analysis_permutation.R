## Genetic correlation permutation test
##
## Use resampling to determine the null distribution of the genetic correlation estimated
## using a bi-variate mixed model and genomic relationship matrix. The stability
## estimates will be permuted, but the mean will not be
##
## Author: Jeff Neyhart
## Last modified: July 20, 2018
##

# Run the source script - MSI
repo_dir <- "/panfs/roc/groups/6/smithkp/neyha001/QTLMapping/S2MET_Mapping/"
source(file.path(repo_dir, "source_MSI.R"))

# repo_dir <- getwd()
# source(file.path(repo_dir, "source.R"))

# Load fw data
load(file.path(result_dir, "pheno_mean_fw_results.RData"))
# Load GWAS results
load(file.path(result_dir, "pheno_fw_mean_gwas_results.RData"))


# Relationship matrix
M <- s2tp_genos_imputed
K <- A.mat(X = M, min.MAF = 0, max.missing = 1)

# Significance threshold
alpha <- 0.05

# Number of cores
n_cores <- detectCores()

# Substitute the estimates of stability using just the TP with estimates using
# both the TP and VP
pheno_mean_fw <- pheno_mean_fw_tpvp %>%
  filter(line_name %in% tp)


# Extract the unique stability coefficients
# Then add the breeding program information
pheno_mean_fw1 <- pheno_mean_fw %>%
  left_join(., subset(entry_list, Class == "S2TP", c(Line, Program)), 
            by = c("line_name" = "Line")) %>%
  rename(program = Program)

## Log transform the non-linear stability estimates
pheno_fw_use <- pheno_mean_fw1 %>%
  group_by(trait) %>% 
  mutate(log_delta = log(delta)) %>%
  # Tidy
  gather(term, estimate, b, delta, log_delta) %>% 
  filter(term != "delta")





# Create a df for modeling
pheno_fw_use_tomodel <- pheno_fw_use %>%
  ungroup() %>%
  distinct(trait, line_name, g, term, estimate)


## Use a permutation test to determine significance
n_perm_iter <- 5000

# Set the seed for reproducibility
set.seed(1024)

# Generate the permutations
pheno_fw_use_tomodel_perm <- pheno_fw_use_tomodel %>%
  group_by(trait, term) %>%
  do(permute(data = ., n = n_perm_iter, g, estimate))

# Assign cores and split
pheno_fw_use_tomodel_perm_split <- pheno_fw_use_tomodel_perm %>%
  ungroup() %>%
  assign_cores(n_cores) %>%
  split(.$core)

# Parallelize
pheno_fw_gen_corr_perm <- mclapply(X = pheno_fw_use_tomodel_perm_split, FUN = function(core_df) {
  # Empty vector
  corr_out <- numeric(nrow(core_df))
  
  # Create the model matrices
  mf <- model.frame(g ~ line_name, as.data.frame(core_df$perm[[1]]))
  Z <- t(model.matrix(~ -1 + line_name, mf)) %>%
    `rownames<-`(., colnames(K))
  X <- t(model.matrix(~ 1, mf))
  
  # Iterate
  for (i in seq_along(corr_out)) {
    Y <- t(as.matrix(as.data.frame(core_df$perm[[i]])[,c("g", "estimate")]))
    fit <- emmremlMultivariate(Y = Y, X = X, Z = Z, K = K)
    vcovG <- fit$Vg
    
    # fit <- sommer::mmer(Y = t(Y), X = t(X), Z = list(gen = list(Z = t(Z), K = K)), silent = TRUE)
    # vcovG <- fit$var.comp$gen
    corr_out[i] <- vcovG[1,2] / prod(sqrt(diag(vcovG)))
    
  }
  
  # Add the vector to the df and return
  core_df %>% 
    mutate(corrG = corr_out) %>% 
    select(trait, term, iter = .id, corrG)
  
}, mc.cores = n_cores)


## Collapse
pheno_fw_gen_corr_perm <- bind_rows(pheno_fw_gen_corr_perm)


## Use the significant GWAS results to correct for large-effect QTL
# Parallelize
pheno_fw_gen_corr_qtl_perm <- mclapply(X = pheno_fw_use_tomodel_perm_split, FUN = function(core_df) {
  # Empty vector
  corr_out <- numeric(nrow(core_df))
  
  # Create the model matrices
  mf <- model.frame(g ~ line_name, as.data.frame(core_df$perm[[1]]))
  Z <- t(model.matrix(~ -1 + line_name, mf)) %>%
    `rownames<-`(., colnames(K))
  Xmu <- model.matrix(~ 1, mf)
  
  
  # Iterate
  for (i in seq_along(corr_out)) {
    
    y <- as.matrix(as.data.frame(core_df$perm[[i]])[,c("g", "estimate")])
    # Y <- t(scale(y))
    Y <- t(y)
    
    ## Get the significant SNPs for this trait/term combination
    snps <- subset(gwas_mlmm_final_tpvp, model == "QG" & trait == core_df$trait[i] & coef %in% c("g", core_df$term[i]), 
                   term, drop = T)
    
    X_snp <- M[,snps,drop = FALSE]
    X <- t(cbind(Xmu, X_snp))
    
    fit <- emmremlMultivariate(Y = Y, X = X, Z = Z, K = K)
    vcovG <- fit$Vg
    
    # fit <- sommer::mmer(Y = t(Y), X = t(X), Z = list(gen = list(Z = t(Z), K = K)), silent = TRUE)
    # vcovG <- fit$var.comp$gen
    corr_out[i] <- vcovG[1,2] / prod(sqrt(diag(vcovG)))
    
  }
  
  # Add the vector to the df and return
  core_df %>% 
    mutate(corrG = corr_out) %>% 
    select(trait, term, iter = .id, corrG)
  
}, mc.cores = n_cores)


## Collapse
pheno_fw_gen_corr_qtl_perm <- bind_rows(pheno_fw_gen_corr_qtl_perm)





# Save
save_file <- file.path(result_dir, "stability_correlation_permutation.RData")
save("pheno_fw_gen_corr_perm", "pheno_fw_gen_corr_qtl_perm", file = save_file)

