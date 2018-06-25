## Genomic prediction of FW sensitivity coefficients
## 
## Conduct genomic prediction and cross-validation of FW sensitivity coefficients
## For CV, use both GBS markers and BOPA markers. Also conduct predictions using the 
## genotypic means
## 
## Author: Jeff Neyhart
## Last updated: June 25, 2018
## 


# Run the source script - MSI
repo_dir <- "/panfs/roc/groups/6/smithkp/neyha001/QTLMapping/S2MET_Mapping/"
source(file.path(repo_dir, "source_MSI.R"))

# # Run the source script - local
# repo_dir <- getwd()
# source(file.path(repo_dir, "source.R"))

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


# Set the number of CV iterations and the number of K folds
n_cv_iter <- 100
n_k <- 7

# Create the cv parts
# Use set.seed for reproducibility
set.seed(1055)

cv_parts <- pheno_mean_fw_tomodel %>%
  group_by(trait, coef) %>%
  do(folds = rerun(.n = n_cv_iter, crossv_kfold(data = ., k = n_k)))

# Create a df
cv_parts_df <- cv_parts %>% 
  unnest() %>%
  group_by(trait, coef) %>%
  mutate(iter = seq(n())) %>% 
  ungroup()

# Split by core
cv_parts_df_core <- cv_parts_df %>% 
  assign_cores(n_core = n_cores) %>% 
  split(.$core)

# Run the cross-validation using the whole marker set
cv_results <- mclapply(X = cv_parts_df_core, FUN = function(core_df) {
  
  cv_out <- map_dbl(core_df$folds, ~as.list(.) %>%
                      pmap(., ~predict_RR(train = ..1, test = ..2, K = K)) %>% 
                      bind_rows() %>%
                      {cor(.$value, .$pred_value)} )
  
  # Add to the cor_df and return
  core_df %>% 
    mutate(acc = cv_out) %>% 
    select(-folds, -core)
  
}, mc.cores = n_cores)

# Bind the list elements together
all_marker_cv_results <- bind_rows(cv_results)




## For the top-rank markers, evenly-spaced markers, and top-rank/evenly spaced, run cross-validation

## Evenly-spaced markers (esm)
# Create relationship matrices
K_esm <- evenly_spaced_markers %>% 
  map("marker") %>% 
  map(~A.mat(X = M[,.], min.MAF = 0, max.missing = 1)) 

# Iterate over these K matrices
esm_cv_results <- structure(vector("list", length(K_esm)), names = names(K_esm))

for (i in seq_along(K_esm)) {
  
  K1 <- K_esm[[i]]
  
  # Run the cross-validation
  cv_results <- mclapply(X = cv_parts_df_core, FUN = function(core_df) {
    
    cv_out <- map_dbl(core_df$folds, ~as.list(.) %>%
                        pmap(., ~predict_RR(train = ..1, test = ..2, K = K1)) %>% 
                        bind_rows() %>%
                        {cor(.$value, .$pred_value)} )
    
    # Add to the cor_df
    core_df %>% 
      mutate(acc = cv_out) %>% 
      select(-folds, -core)
    
  }, mc.cores = n_cores)
  
  
  # Bind and add to the list
  esm_cv_results[[i]] <- bind_rows(cv_results)
  
}


## Ranked markers (top)
# Create relationship matrices
K_top_list <- top_rank_markers$marker_subsets %>% 
  map(~map(., "marker") %>%
        map(~A.mat(X = M[,.], min.MAF = 0, max.missing = 1)) )

## Add K matrices and data to the df, then split on cores
cv_parts_df_core <- top_rank_markers %>% 
  mutate(K_mat = K_top_list) %>%
  left_join(cv_parts_df, .) %>% 
  select(-marker_subsets) %>%   
  assign_cores(n_core = n_cores) %>% 
  split(.$core)

# Run the cross-validation
cv_results <- mclapply(X = cv_parts_df_core, FUN = function(core_df) {
  
  # Create an empty vector for storing the accuracy results for each marker subset group
  acc_out <- vector("list", nrow(core_df))
  
  # Seq along the rows
  for (i in seq(nrow(core_df))) {
    
    # Extract the K matrices and the folds
    folds <- core_df$folds[[i]]
    K_mats <- core_df$K_mat[[i]]
    
    # Map over the Ks
    K_acc_out <- map_dbl(K_mats, function(K) pmap_df(folds, ~predict_RR(train = ..1, test = ..2, K = K)) %>%
                       {cor(.$value, .$pred_value)} )
    
    ## Create a df and add to the list
    acc_out[[i]] <- K_acc_out %>% 
      data.frame(nmar = names(.), acc = ., row.names = NULL, stringsAsFactors = FALSE)
      
  }
  
  ## Add the accuracy list to the core_df and return
  core_df %>% 
    mutate(acc = acc_out) %>% 
    select(-folds, -K_mat, -core)
  
}, mc.cores = n_cores)

# Bind the list elements together
top_marker_cv_results <- bind_rows(cv_results)

  

## Ranked, evenly-spaced markers (tesm)
# Create relationship matrices
K_tesm_list <- top_rank_evenly_spaced_markers$marker_subsets %>% 
  map(~map(., "marker") %>%
        map(~A.mat(X = M[,.], min.MAF = 0, max.missing = 1)) )

## Add K matrices and data to the df, then split on cores
cv_parts_df_core <- top_rank_evenly_spaced_markers %>% 
  mutate(K_mat = K_tesm_list) %>%
  left_join(cv_parts_df, .) %>% 
  select(-marker_subsets) %>%   
  assign_cores(n_core = n_cores) %>% 
  split(.$core)

# Run the cross-validation
cv_results <- mclapply(X = cv_parts_df_core, FUN = function(core_df) {
  
  # Create an empty vector for storing the accuracy results for each marker subset group
  acc_out <- vector("list", nrow(core_df))
  
  # Seq along the rows
  for (i in seq(nrow(core_df))) {
    
    # Extract the K matrices and the folds
    folds <- core_df$folds[[i]]
    K_mats <- core_df$K_mat[[i]]
    
    # Map over the Ks
    K_acc_out <- map_dbl(K_mats, function(K) pmap_df(folds, ~predict_RR(train = ..1, test = ..2, K = K)) %>%
    {cor(.$value, .$pred_value)} )
    
    ## Create a df and add to the list
    acc_out[[i]] <- K_acc_out %>% 
      data.frame(nmar = names(.), acc = ., row.names = NULL, stringsAsFactors = FALSE)
    
  }
  
  ## Add the accuracy list to the core_df and return
  core_df %>% 
    mutate(acc = acc_out) %>% 
    select(-folds, -K_mat, -core)
  
}, mc.cores = n_cores)

# Bind the list elements together
tesm_marker_cv_results <- bind_rows(cv_results)


## Top plastic markers (plas)
# Only use the top plastic markers (not the top intercept markers)
top_plastic_markers1 <- top_plastic_markers %>% 
  filter(coef == "c") %>%
  select(-coef)


# Create relationship matrices
K_plas_list <- top_plastic_markers1$marker_subsets %>% 
  map(~map(., "marker") %>%
        map(~A.mat(X = M[,.], min.MAF = 0, max.missing = 1)) )

## Add K matrices and data to the df, then split on cores
cv_parts_df_core <- top_plastic_markers1 %>% 
  mutate(K_mat = K_plas_list) %>%
  left_join(cv_parts_df, .) %>% 
  select(-marker_subsets) %>%   
  assign_cores(n_core = n_cores) %>% 
  split(.$core)

# Run the cross-validation
cv_results <- mclapply(X = cv_parts_df_core, FUN = function(core_df) {
  
  # Create an empty vector for storing the accuracy results for each marker subset group
  acc_out <- vector("list", nrow(core_df))
  
  # Seq along the rows
  for (i in seq(nrow(core_df))) {
    
    # Extract the K matrices and the folds
    folds <- core_df$folds[[i]]
    K_mats <- core_df$K_mat[[i]]
    
    # Map over the Ks
    K_acc_out <- map_dbl(K_mats, function(K) pmap_df(folds, ~predict_RR(train = ..1, test = ..2, K = K)) %>%
    {cor(.$value, .$pred_value)} )
    
    ## Create a df and add to the list
    acc_out[[i]] <- K_acc_out %>% 
      data.frame(nmar = names(.), acc = ., row.names = NULL, stringsAsFactors = FALSE)
    
  }
  
  ## Add the accuracy list to the core_df and return
  core_df %>% 
    mutate(acc = acc_out) %>% 
    select(-folds, -K_mat, -core)
  
}, mc.cores = n_cores)

# Bind the list elements together
plas_marker_cv_results <- bind_rows(cv_results)




## Random marker subsets

## We have a list of randomly sampled markers for each marker subset size
## Generate 10 CV iteration
n_cv_iter <- 10
n_samples <- length(random_markers[[1]])

# Create the relationship matrices
K_rand_list <- random_markers %>% 
  map(~map(., "marker") %>% map(~A.mat(X = M[,.], min.MAF = 0, max.missing = 1)) )

# Create a df
K_rand_df <- K_rand_list %>% 
  as_data_frame() %>% 
  mutate(iter = seq(nrow(.))) %>% 
  gather(nmar, K_mat, -iter)


## Create CV samples
# Use set.seed for reproducibility
set.seed(1138)

cv_parts <- pheno_mean_fw_tomodel %>%
  group_by(trait, coef) %>%
  nest(line_name, value) %>%
  mutate(folds_list = map(data, ~replicate(n = n_samples, replicate(n = n_cv_iter, crossv_kfold(data = ., k = n_k), simplify = FALSE), 
                                           simplify = FALSE)))
  
# Create a df
cv_parts_df <- cv_parts %>% 
  select(-data) %>% 
  unnest() %>% 
  group_by(trait, coef) %>%
  mutate(iter = seq(n())) %>% 
  ungroup()


# Split by core
cv_parts_df_core <- cv_parts_df %>% 
  assign_cores(n_core = n_cores) %>% 
  split(.$core)


# Run the cross-validation
cv_results <- mclapply(X = cv_parts_df_core, FUN = function(core_df) {
  
  ## Add the relationship matrix samples
  core_df1 <- left_join(core_df, K_rand_df)
  
  # Create an empty vector for storing the accuracy results for each marker subset group
  acc_out <- vector("numeric", nrow(core_df1))
  
  # Seq along the rows
  for (i in seq(nrow(core_df1))) {
    
    # Extract the K matrices and the folds
    folds_list <- core_df1$folds_list[[i]]
    K_mat <- core_df1$K_mat[[i]]
    
    # Map over the folds
    fold_acc_out <- map_dbl(folds_list, ~pmap_df(., ~predict_RR(train = ..1, test = ..2, K = K_mat)) %>%
                           {cor(.$value, .$pred_value)} )
    
    ## Take the mean and add it to the list
    acc_out[[i]] <- mean(fold_acc_out)
    
  }
  
  ## Add the accuracy list to the core_df and return
  core_df1 %>% 
    mutate(acc = acc_out) %>% 
    select(-folds_list, -K_mat, -core)
  
}, mc.cores = n_cores)

# Bind the list elements together
rand_marker_cv_results <- bind_rows(cv_results)


## Save everything
save_file <- file.path(result_dir, "stability_mean_crossv_results.RData")
save("all_marker_cv_results", "esm_cv_results", "top_marker_cv_results", "tesm_marker_cv_results", 
     "plas_marker_cv_results", "rand_marker_cv_results", file = save_file)

