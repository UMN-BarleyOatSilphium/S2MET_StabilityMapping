## Variance explained by different marker subsets and groups
## 
## Calculate the variance explained by different marker groups and numbers of markers
## 
## Author: Jeff Neyhart
## Last updated: June 25, 2018
## 


# Run the source script - local
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))
library(modelr)
library(parallel)

# Load the FW results
load(file.path(result_dir, "pheno_mean_fw_results.RData"))
# Load the marker subsets
load(file.path(result_dir, "marker_subsets.RData"))

# Cores
n_cores <- detectCores()


# Rename the marker matrix
M <- S2TP_imputed_multi_genos_mat[tp_geno,]
# Overall K
K <- A.mat(X = M, min.MAF = 0, max.missing = 1)


# Create a tidy dataset to model
pheno_mean_fw_tomodel <- pheno_mean_fw %>% 
  distinct(trait, line_name, g, b, delta) %>%
  mutate(line_name = as.factor(line_name),
         log_delta = log(delta)) %>%
  select(-delta) %>%
  gather(coef, value, g:log_delta)

# For each trait and coef, create model matrices
pheno_mean_fw_tomodel1 <- pheno_mean_fw_tomodel %>%
  group_by(trait, coef) %>% 
  do({
    mf <- model.frame(value ~ line_name, .)
    y <- model.response(mf)
    X <- model.matrix(~ 1, mf)
    Z <- model.matrix(~ -1 + line_name, mf)
    data_frame(y = list(y), X = list(X), Z = list(Z))
  })
    


## Variance explained by all markers
all_marker_varcomp <- pheno_mean_fw_tomodel1 %>%
  group_by(trait, coef) %>% 
  do({
    df <- .
    
    # Fit the model
    fit <- mixed.solve(y = df$y[[1]], Z = df$Z[[1]], K = K, X = df$X[[1]])
    data.frame(varG = fit$Vu, varR = fit$Ve, loglik = fit$LL)
    
  })


## Evenly-spaced markers (esm)
# Create relationship matrices
K_esm <- evenly_spaced_markers %>% 
  map("marker") %>% 
  map(~A.mat(X = M[,.], min.MAF = 0, max.missing = 1)) 

## Variance explained
esm_marker_varcomp <- pheno_mean_fw_tomodel1 %>%
  group_by(trait, coef) %>% 
  do({
    df <- .
    
    # Map over the K matrices
    out <- map_df(K_esm, ~{
      # Fit the model
      fit <- mixed.solve(y = df$y[[1]], Z = df$Z[[1]], K = ., X = df$X[[1]])
      data.frame(varG = fit$Vu, varR = fit$Ve, loglik = fit$LL)
    })
    
    # Add marker numbers
    mutate(out, nmar = names(K_esm))
    
  })




## Ranked markers (top)
# Create relationship matrices
K_top_list <- top_rank_markers %>%
  map(~group_by(., trait, coef) %>%
        do(K_mat = A.mat(X = M[,.$marker], min.MAF = 0, max.missing = 1)) %>%
        ungroup() )


## Variance explained
top_marker_varcomp <- pheno_mean_fw_tomodel1 %>%
  group_by(trait, coef) %>% 
  do({
    df <- .
    
    # Map over the K matrices
    out <- map(K_top_list, ~subset(., trait == df$trait & coef == df$coef, K_mat, drop = TRUE)[[1]]) %>%
      map_df(~{
        # Fit the model
        fit <- mixed.solve(y = df$y[[1]], Z = df$Z[[1]], K = ., X = df$X[[1]])
        data.frame(varG = fit$Vu, varR = fit$Ve, loglik = fit$LL)
      })
    
    # Add marker numbers
    mutate(out, nmar = names(K_top_list))
    
  }) %>% ungroup()


## Ranked, evenly-spaced markers (tesm)
# Create relationship matrices
K_tesm_list <- top_rank_evenly_spaced_markers %>%
  map(~group_by(., trait, coef) %>%
        do(K_mat = A.mat(X = M[,.$marker], min.MAF = 0, max.missing = 1)) %>%
        ungroup() )


## Variance explained
tesm_marker_varcomp <- pheno_mean_fw_tomodel1 %>%
  group_by(trait, coef) %>% 
  do({
    df <- .
    
    # Map over the K matrices
    out <- map(K_tesm_list, ~subset(., trait == df$trait & coef == df$coef, K_mat, drop = TRUE)[[1]]) %>%
      map_df(~{
        # Fit the model
        fit <- mixed.solve(y = df$y[[1]], Z = df$Z[[1]], K = ., X = df$X[[1]])
        data.frame(varG = fit$Vu, varR = fit$Ve, loglik = fit$LL)
      })
    
    # Add marker numbers
    mutate(out, nmar = names(K_tesm_list))
    
  }) %>% ungroup()




## Top plastic markers (plas)
# Only use the top plastic markers (not the top intercept markers)
top_plastic_markers1 <- top_plastic_markers %>% 
  filter(coef == "c") %>%
  select(-coef)


# Create relationship matrices
K_plas_list <- top_plastic_markers1$marker_subsets %>% 
  map(~map(., "marker") %>%
        map(~A.mat(X = M[,.], min.MAF = 0, max.missing = 1)) )


# Create a df
K_plas_df <- top_plastic_markers1 %>% 
  mutate(K_mat = K_plas_list)

## Variance explained
plas_marker_varcomp <- pheno_mean_fw_tomodel1 %>%
  left_join(., K_plas_df) %>%
  group_by(trait, coef) %>% 
  do({
    df <- .
    
    # Map over the K matrices
    out <- map_df(df$K_mat[[1]], ~{
      # Fit the model
      fit <- mixed.solve(y = df$y[[1]], Z = df$Z[[1]], K = ., X = df$X[[1]])
      data.frame(varG = fit$Vu, varR = fit$Ve, loglik = fit$LL)
    })
    
    # Add marker numbers
    mutate(out, nmar = names(K_esm))
    
  })


## Random markers (rand)

# Create the relationship matrices
K_rand_list <- random_markers %>% 
  # map(~.[seq(n_samples)]) %>%
  map(~map(., "marker") %>% map(~A.mat(X = M[,.], min.MAF = 0, max.missing = 1)) )

# Create a df
K_rand_df <- K_rand_list %>% 
  as_data_frame() %>% 
  mutate(iter = seq(nrow(.))) %>% 
  gather(nmar, K_mat, -iter)

## Add cores and split
K_rand_df_split <- K_rand_df %>%
  assign_cores(n_cores) %>%
  split(.$core)


# Parallelize the model fittings
rand_marker_varcomp_out <- mclapply(X = K_rand_df_split, FUN = function(core_df) {
  
  pheno_mean_fw_tomodel1 %>%
    group_by(trait, coef) %>% 
    do({
      df <- .
      
      # Map over the K matrices
      out <- map(core_df$K_mat, ~{
        # Fit the model
        fit <- mixed.solve(y = df$y[[1]], Z = df$Z[[1]], K = ., X = df$X[[1]])
        data.frame(varG = fit$Vu, varR = fit$Ve, loglik = fit$LL)
      })
      
      core_df %>% 
        mutate(out = out) %>% 
        select(-K_mat, -core) %>%
        unnest()
      
    })
  
}, mc.cores = n_cores)


# Bind the list elements together
rand_marker_varcomp <- bind_rows(rand_marker_varcomp_out) %>%
  ungroup()



## Save everything
save_file <- file.path(result_dir, "stability_mean_marker_varcomp_results.RData")
save("all_marker_varcomp", "esm_marker_varcomp", "top_marker_varcomp", "tesm_marker_varcomp", 
     "plas_marker_varcomp", "rand_marker_varcomp", file = save_file)

