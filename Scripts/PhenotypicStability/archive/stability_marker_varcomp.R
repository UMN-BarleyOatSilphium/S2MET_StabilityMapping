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
load(file.path(result_dir, "marker_mean_fw_results.RData"))

# Cores
n_cores <- detectCores()


# Rename the marker matrix
M <- s2tp_genos_imputed[tp_geno,]
# Overall K
K <- A.mat(X = M, min.MAF = 0, max.missing = 1)

# Use the TP-VP results
pheno_mean_fw <- pheno_mean_fw_tpvp %>%
  filter(line_name %in% tp_geno)


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
    
  }) %>% ungroup()







top_plastic_markers <- marker_mean_fw_ann %>% 
  filter(annotation == "plastic") %>%
  split(.$trait) %>%
  map("marker")

## Top plastic markers (plas)

# Create relationship matrices
K_plas_list <- top_plastic_markers %>%
  map(~A.mat(X = M[,.], min.MAF = 0, max.missing = 1))

# Create a df
K_plas_df <- K_plas_list %>%
  data_frame(trait = names(.), K_mat = .)


## Create samples of stable markers
stable_markers <- marker_mean_fw_ann %>% 
  filter(annotation == "stable") %>%
  split(.$trait) %>%
  map("marker")

n_samples <- 100

stable_markers_samples <- map2(.x = top_plastic_markers, .y = stable_markers, ~{
  replicate(n = n_samples, expr = sample(x = .y, size = length(.x)), simplify = FALSE) })

## Create relationship matrices
K_stab_list <- stable_markers_samples %>%
  map(., ~map(., ~A.mat(X = M[,.], min.MAF = 0, max.missing = 1)))

# Create a df
K_stab_df <- K_stab_list %>%
  data_frame(trait = names(.), K_mat = .)


## Variance explained
plas_stab_marker_varcomp <- pheno_mean_fw_tomodel1 %>%
  left_join(., K_plas_df) %>%
  left_join(., K_stab_df, by = c("trait")) %>%
  group_by(trait, coef) %>% 
  do({
    df <- .
    
    # Iterate over the list of stab marker matrices
    fit_list <- df$K_mat.y[[1]] %>%
      map(~{
        K2 <- .
        
        EMMREML::emmremlMultiKernel(y = df$y[[1]], X = df$X[[1]], Zlist = list(df$Z[[1]], df$Z[[1]]),
                                    Klist = list(df$K_mat.x[[1]], K2))
        
      })
    
    # Organize
    fit_list %>% 
      map_df(~{data.frame(plastic = .$Vu * .$weights[1], stable = .$Vu * .$weights[2], residual = .$Ve)}) %>%
      mutate(iter = seq(nrow(.))) %>%
      select(iter, names(.))

  })





## Save everything
save_file <- file.path(result_dir, "stability_mean_marker_varcomp_results.RData")
save("all_marker_varcomp", "plas_stab_marker_varcomp", file = save_file)

