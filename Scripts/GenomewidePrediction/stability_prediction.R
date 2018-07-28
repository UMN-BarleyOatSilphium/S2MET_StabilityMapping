## Genomic prediction of FW sensitivity coefficients
## 
## Conduct genomic prediction and cross-validation of FW sensitivity coefficients
## For CV, use both GBS markers and BOPA markers. Also conduct predictions using the 
## genotypic means
## 
## Author: Jeff Neyhart
## Last updated: June 26, 2018
## 


# Run the source script - local
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

# Other packages
library(modelr)
library(parallel)
library(lme4)
library(broom)
library(cowplot)

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

# Colors for the coefficients
colors_use <- set_names(umn_palette(2)[3:5], coef_replace)


# Number of cv iterations
n_cv_iter <- 100


## TP cross-validation
p_train <- 0.60

# Use the TP + VP data
pheno_mean_fw <- pheno_mean_fw_tpvp %>%
  filter(line_name %in% tp)


# Create a tidy data.frame for modeling
pheno_mean_fw_tomodel <- pheno_mean_fw %>%
  distinct(trait, line_name, g, b, delta) %>%
  mutate(log_delta = log(delta)) %>%
  select(-delta) %>%
  gather(coef, value, g:log_delta)


## Create the cross-validation samples
pheno_mean_fw_tomodel_cv <- pheno_mean_fw_tomodel %>% 
  mutate(line_name = as.factor(line_name)) %>%
  group_by(trait, coef) %>% 
  do(cv_sample = crossv_mc(data = ., n = n_cv_iter, test = 1 - p_train)) %>%
  ungroup()

## Run cross-val
cv_results_alldata <- pheno_mean_fw_tomodel_cv %>%
  group_by(trait, coef) %>%
  do({
    
    df <- .
    
    # Map, predict, and measure accuracy
    acc_out <- df$cv_sample[[1]] %>% 
      pmap(~predict_RR(train = ..1, test = ..2, K = K)) %>%
      map_dbl(~cor(.$value, .$pred_value))
    
    df %>% 
      unnest() %>%
      mutate(iter = seq(n()), acc = acc_out) %>% 
      select(-train:-.id)
    
  })


  



## Run similar 60:40 CV, but instead predict the value of the unobserved genotypes in
## each environment
pheno_mean_fw_tomodel <- pheno_mean_fw %>%
  select(trait, line_name, environment, value, std_error, g, b, h)


## Create the cross-validation samples
pheno_mean_fw_tomodel_cv <- pheno_mean_fw_tomodel %>% 
  mutate(line_name = as.factor(line_name)) %>%
  group_by(trait) %>% 
  do(cv_sample = {
    df <- . 
    
    # Create a set of training lines (we assume we know the mean and regression slope)
    sample_train_line_name <- replicate(n = n_cv_iter, sample_frac(tbl = distinct(df, line_name), size = p_train), simplify = FALSE)
    # Create a set of training environments
    sample_train_environment <- replicate(n = n_cv_iter, sample_frac(tbl = distinct(df, environment), size = p_train), simplify = FALSE)
    
    # Create training sets of the mean and stability of the TP
    training_data <- list(sample_train_line_name, sample_train_environment) %>%
      pmap(~filter(df, line_name %in% .x[[1]], environment %in% .y[[1]]))
    
    
    # Create testing sets of the environment-specific observations of the testing pop
    testing_data <- training_data %>%
      map(~filter(df, !line_name %in% unique(.$line_name), !environment %in% unique(.$environment)))
    
    # Return the training and testing sets
    data_frame(train = training_data, test = testing_data)
    
  })
    

## Extract the line names for creating factors
line_name_levels <- unique(pheno_mean_fw_tomodel$line_name)

## Run cross-val
cv_results_env <- pheno_mean_fw_tomodel_cv %>%
  group_by(trait) %>%
  do({
    df <- .
    
    ## Iterate over the training data and calculate genotype means and stabilities
    training_data_topred <- df %>% 
      unnest() %>% 
      pull(train) %>% 
      map(~{
        # Calculate the genotype and environmental mean
        gh <-select(., -g:-h) %>% 
          calc_gh()
        # Calculate stability
        stab <- gh %>% 
          group_by(line_name) %>% 
          do(calc_stability(., remove.outliers = FALSE)) %>% 
          filter(type == "outliers") %>% 
          select(line_name, b, delta)
        
        # Join the data
        gh %>% 
          left_join(stab, by = "line_name") %>%
          distinct(trait, line_name, g, b, delta) %>%
          mutate(log_delta = log(delta)) %>%
          select(-delta)
        
      })
    
    # Do the same for the testing set
    testing_data_topred <- df %>% 
      unnest() %>% 
      pull(test) %>% 
      map(~{
        # Calculate the genotype and environmental mean
        gh <- select(., -g:-h) %>% 
          calc_gh()
        # Calculate stability
        stab <- gh %>% 
          group_by(line_name) %>% 
          do(calc_stability(., remove.outliers = FALSE)) %>% 
          filter(type == "outliers") %>% 
          select(line_name, b, delta)
        
        # Join the data
        gh %>% 
          left_join(stab, by = "line_name") %>%
          mutate(log_delta = log(delta)) %>%
          select(-delta)
        
      })
    
    # Predict and assess accuracy
    acc_out <- list(training_data_topred, testing_data_topred) %>%
      pmap_df(~{
        train <- mutate(.x, line_name = factor(line_name, levels = line_name_levels))
       
        g_pred_out <- mixed.solve(y = train$g, model.matrix(~ -1 + line_name, train), K = K)$u %>%
          data_frame(line_name = names(.), pred_g = .)
        b_pred_out <- mixed.solve(y = train$b, model.matrix(~ -1 + line_name, train), K = K)$u %>%
          data_frame(line_name = names(.), pred_b = .)
        d_pred_out <- mixed.solve(y = train$log_delta, model.matrix(~ -1 + line_name, train), K = K)$u %>%
          data_frame(line_name = names(.), pred_log_delta = .)
        
        # Predict the mean and stability of the testing set
        pred_gb <- reduce(list(.y, g_pred_out, b_pred_out, d_pred_out), left_join, by = "line_name") %>%
          distinct(line_name, g, b, log_delta, pred_g, pred_b, pred_log_delta) %>% 
          summarize(acc_pred_g = cor(g, pred_g), 
                    acc_pred_b = cor(b, pred_b),
                    acc_pred_log_delta = cor(log_delta, pred_log_delta))
        
        # Predict the phenotype in each environment
        pred_pheno <- reduce(list(.y, g_pred_out, b_pred_out), left_join, by = "line_name") %>%
          mutate(pred_value_gb = pred_g + (pred_b * h),
                 pred_value_g = pred_g) %>% # Also just use the genotype mean 
          group_by(environment) %>% 
          summarize(acc_gb = cor(value, pred_value_gb),
                    acc_g = cor(value, pred_value_g))
        
        data_frame(pred_gb = list(pred_gb), pred_pheno = list(pred_pheno))
      })
    
    
    # Add the accuracy results to the df
    df %>% 
      unnest() %>% 
      bind_cols(., acc_out) %>%
      mutate(iter = seq(n())) %>% 
      select(trait, iter, pred_gb, pred_pheno)
    
  })




#### Predictions of the VP


# Rename the marker matrix
M <- s2_imputed_mat[c(tp_geno1, vp_geno1),]
# Overall K
K <- A.mat(X = M, min.MAF = 0, max.missing = 1)


# Create a tidy dataset to model
pheno_mean_fw_tomodel <- pheno_mean_fw_tpvp %>% 
  filter(line_name %in% c(tp_geno1, vp_geno1)) %>%
  distinct(trait, line_name, g, b, delta) %>%
  mutate(line_name = as.factor(line_name),
         log_delta = log(delta)) %>%
  select(-delta) %>%
  gather(coef, value, g:log_delta)


## For each trait, predict the mean, slope, and MSE of the validation set
vp_prediction_all_data <- pheno_mean_fw_tomodel %>% 
  mutate(line_name = as.factor(line_name)) %>%
  group_by(trait, coef) %>%
  do({
    df <- .
    
    # Separate the training from validation
    df_train <- filter(df, line_name %in% tp_geno1)
    mf <- model.frame(value ~ line_name, df_train)
    y <- model.response(mf)
    Z <- model.matrix(~ -1 + line_name, mf)
    
    # Fit the model
    fit <- mixed.solve(y = y, Z = Z, K = K)
    
    # Return the BLUPs
    blups <- data.frame(line_name = names(fit$u), pred_value = fit$u, row.names = NULL, stringsAsFactors = FALSE)
    # Just the vp
    filter(df, line_name %in% vp_geno1) %>%
      left_join(., blups, by = "line_name")
    
  })




## Now sample environments for the training set and predict values of the VP
pheno_mean_fw_tomodel <- pheno_mean_fw_tpvp %>%
  filter(line_name %in% c(tp_geno1, vp_geno1)) %>%
  select(trait, line_name, environment, value, std_error, g, b, h)

# Create vectors of environments with at least one training individual or at least one validation
# individual
train_environments <- pheno_mean_fw_tomodel %>% 
  group_by(trait, environment) %>% 
  filter(sum(line_name %in% tp_geno1) >= 1) %>% 
  distinct(trait, environment) %>%
  ungroup()

val_environments <- pheno_mean_fw_tomodel %>% 
  group_by(trait, environment) %>% 
  filter(sum(line_name %in% vp_geno1) >= 1) %>% 
  distinct(trait, environment) %>%
  ungroup()


## Create the parent-offspring validation (pov) samples
pheno_mean_fw_tomodel_pov <- pheno_mean_fw_tomodel %>% 
  mutate(line_name = as.factor(line_name)) %>%
  group_by(trait) %>% 
  do(cv_sample = {
    df <- . 
    
    # Create a set of training environments
    sample_train_environment <- replicate(n = n_cv_iter, sample_frac(tbl = filter(train_environments, trait == unique(df$trait)), 
                                                                     size = p_train), simplify = FALSE)
    
    # Sample the remaining validation environments
    sample_val_environment <- sample_train_environment %>%
      map(~filter(val_environments, trait == unique(df$trait), !environment %in% .$environment))
    
    
    # Create training sets of the mean and stability of the TP
    training_data <- sample_train_environment %>%
      map(~filter(df, line_name %in% tp_geno1, environment %in% .$environment))
    
    
    # Create testing sets of the environment-specific observations of the VP
    validation_data <- sample_val_environment %>%
      map(~filter(df, line_name %in% vp_geno1, environment %in% .$environment))
    
    # Return the training and testing sets
    data_frame(train = training_data, test = validation_data)
    
  })


## Extract the line names for creating factors
line_name_levels <- unique(pheno_mean_fw_tomodel$line_name)

## Run parent-offpring validation
vp_prediction_env <- pheno_mean_fw_tomodel_pov %>%
  group_by(trait) %>%
  do({
    df <- .
    
    ## Iterate over the training data and calculate genotype means and stabilities
    training_data_topred <- df %>% 
      unnest() %>% 
      pull(train) %>% 
      map(~{
        # Calculate the genotype and environmental mean
        gh <-select(., -g:-h) %>% 
          calc_gh()
        # Calculate stability
        stab <- gh %>% 
          group_by(line_name) %>% 
          do(calc_stability(., remove.outliers = FALSE)) %>% 
          filter(type == "outliers") %>% 
          select(line_name, b, delta)
        
        # Join the data
        gh %>% 
          left_join(stab, by = "line_name") %>%
          distinct(trait, line_name, g, b)
        
      })
    
    # Do the same for the testing set
    testing_data_topred <- df %>% 
      unnest() %>% 
      pull(test) %>% 
      map(~{
        # Calculate the genotype and environmental mean
        gh <- select(., -g:-h) %>% 
          calc_gh()
        # Calculate stability
        stab <- gh %>% 
          group_by(line_name) %>% 
          do(calc_stability(., remove.outliers = FALSE)) %>% 
          filter(type == "outliers") %>% 
          select(line_name, b, delta)
        
        # Join the data
        gh %>% 
          left_join(stab, by = "line_name")
        
      })
    
    # Predict and assess accuracy
    acc_out <- list(training_data_topred, testing_data_topred) %>%
      pmap_df(~{
        train <- mutate(.x, line_name = factor(line_name, levels = line_name_levels))
        
        g_pred_out <- mixed.solve(y = train$g, model.matrix(~ -1 + line_name, train), K = K)$u %>%
          data_frame(line_name = names(.), pred_g = .)
        b_pred_out <- mixed.solve(y = train$b, model.matrix(~ -1 + line_name, train), K = K)$u %>%
          data_frame(line_name = names(.), pred_b = .)
        
        # Predict the mean and stability of the testing set
        pred_gb <- reduce(list(.y, g_pred_out, b_pred_out), left_join, by = "line_name") %>%
          distinct(line_name, g, b, pred_g, pred_b) %>% 
          summarize(acc_pred_g = cor(g, pred_g), 
                    acc_pred_b = cor(b, pred_b))
        
        # Predict the phenotype in each environment
        pred_pheno <- reduce(list(.y, g_pred_out, b_pred_out), left_join, by = "line_name") %>%
          mutate(pred_value_gb = pred_g + (pred_b * h),
                 pred_value_g = pred_g) %>% # Also just use the genotype mean 
          group_by(environment) %>% 
          summarize(acc_gb = cor(value, pred_value_gb),
                    acc_g = cor(value, pred_value_g))
        
        data_frame(pred_gb = list(pred_gb), pred_pheno = list(pred_pheno))
      })
    
    
    # Add the accuracy results to the df
    df %>% 
      unnest() %>% 
      bind_cols(., acc_out) %>%
      mutate(iter = seq(n())) %>% 
      select(trait, iter, pred_gb, pred_pheno)
    
  })



## Save this
save_file <- file.path(result_dir, "genomewide_prediction_results.RData")
save("cv_results_alldata", "cv_results_env", "vp_prediction_all_data", "vp_prediction_env",
     file = save_file)




## Use marker subsets to predict stability  
# Load the marker subsets
load(file.path(result_dir, "marker_subsets_tpvp.RData"))

## Top markers
top_rank_markers_predictions <- top_rank_markers_tpvp %>%
  map(function(marker_df) {
    
    pheno_mean_fw_tomodel %>% 
      mutate(line_name = as.factor(line_name)) %>%
      group_by(trait, coef) %>%
      do({
        df <- .
        
        # Grab the markers
        markers_use <- subset(marker_df, trait == unique(df$trait) & coef == unique(df$coef), marker, drop = TRUE)
        # Create relationship matrix
        K_mat <- A.mat(X = M[,markers_use], min.MAF = 0, max.missing = 1)
        
        # Separate the training from validation
        df_train <- filter(df, line_name %in% tp_geno1)
        mf <- model.frame(value ~ line_name, df_train)
        y <- model.response(mf)
        Z <- model.matrix(~ -1 + line_name, mf)
        
        # Fit the model
        fit <- mixed.solve(y = y, Z = Z, K = K_mat)
        
        # Return the BLUPs
        blups <- data.frame(line_name = names(fit$u), pred_value = fit$u, row.names = NULL, stringsAsFactors = FALSE)
        # Just the vp
        filter(df, line_name %in% vp_geno1) %>%
          left_join(., blups, by = "line_name")
        
      })
    
  })
    


## Top, evenly-spaced markers
tesm_markers_predictions <- top_rank_evenly_spaced_markers_tpvp %>%
  map(function(marker_df) {
    
    pheno_mean_fw_tomodel %>% 
      mutate(line_name = as.factor(line_name)) %>%
      group_by(trait, coef) %>%
      do({
        df <- .
        
        # Grab the markers
        markers_use <- subset(marker_df, trait == unique(df$trait) & coef == unique(df$coef), marker, drop = TRUE)
        # Create relationship matrix
        K_mat <- A.mat(X = M[,markers_use], min.MAF = 0, max.missing = 1)
        
        # Separate the training from validation
        df_train <- filter(df, line_name %in% tp_geno1)
        mf <- model.frame(value ~ line_name, df_train)
        y <- model.response(mf)
        Z <- model.matrix(~ -1 + line_name, mf)
        
        # Fit the model
        fit <- mixed.solve(y = y, Z = Z, K = K_mat)
        
        # Return the BLUPs
        blups <- data.frame(line_name = names(fit$u), pred_value = fit$u, row.names = NULL, stringsAsFactors = FALSE)
        # Just the vp
        filter(df, line_name %in% vp_geno1) %>%
          left_join(., blups, by = "line_name")
        
      })
    
  })


## Evenly-spaced markers
esm_markers_predictions <- evenly_spaced_markers_tpvp %>%
  map(function(marker_df) {
    
    # Grab the markers
    markers_use <- subset(marker_df, , marker, drop = TRUE)
    # Create relationship matrix
    K_mat <- A.mat(X = M[,markers_use], min.MAF = 0, max.missing = 1)
    
    pheno_mean_fw_tomodel %>% 
      mutate(line_name = as.factor(line_name)) %>%
      group_by(trait, coef) %>%
      do({
        df <- .
        
        # Separate the training from validation
        df_train <- filter(df, line_name %in% tp_geno1)
        mf <- model.frame(value ~ line_name, df_train)
        y <- model.response(mf)
        Z <- model.matrix(~ -1 + line_name, mf)
        
        # Fit the model
        fit <- mixed.solve(y = y, Z = Z, K = K_mat)
        
        # Return the BLUPs
        blups <- data.frame(line_name = names(fit$u), pred_value = fit$u, row.names = NULL, stringsAsFactors = FALSE)
        # Just the vp
        filter(df, line_name %in% vp_geno1) %>%
          left_join(., blups, by = "line_name")
        
      })
    
  })


## Random markers
# Create a df to parallelize
random_markers_split <- data_frame(nmar = names(random_markers_tpvp), random_markers_tpvp) %>% 
  unnest() %>% 
  group_by(nmar) %>%
  mutate(iter = seq(n())) %>% 
  ungroup() %>%
  assign_cores(n_cores) %>%
  split(.$core)

# Parallelize
rand_marker_predction_out <- mclapply(X = random_markers_split, FUN = function(core_df) {
  
  # Map over the markers
  out <- core_df$random_markers_tpvp %>%
    map("marker") %>%
    map(~{
      
      K_mat <- A.mat(X = M[,.], min.MAF = 0, max.missing = 1)
      
      pheno_mean_fw_tomodel %>% 
        mutate(line_name = as.factor(line_name)) %>%
        group_by(trait, coef) %>%
        do({
          df <- .
          
          # Separate the training from validation
          df_train <- filter(df, line_name %in% tp_geno1)
          mf <- model.frame(value ~ line_name, df_train)
          y <- model.response(mf)
          Z <- model.matrix(~ -1 + line_name, mf)
          
          # Fit the model
          fit <- mixed.solve(y = y, Z = Z, K = K_mat)
          
          # Return the BLUPs
          blups <- data.frame(line_name = names(fit$u), pred_value = fit$u, row.names = NULL, stringsAsFactors = FALSE)
          # Just the vp
          filter(df, line_name %in% vp_geno1) %>%
            left_join(., blups, by = "line_name")
          
        })
      
    })
  
  core_df %>%
    mutate(results = out) %>%
    select(-random_markers_tpvp, -core)
  
}, mc.cores = n_cores)


rand_marker_predctions <- bind_rows(rand_marker_predction_out)


## Save everything
all_marker_prediction_results <- list(
  vp_prediction_all_data = vp_prediction_all_data,
  vp_prediction_loyo = vp_prediction_loyo,
  vp_prediction_rand_env = vp_prediction_rand_env
)

marker_subset_prediction_results <- list(
  top_rank_markers_predictions = top_rank_markers_predictions,
  tesm_markers_predictions = tesm_markers_predictions,
  esm_markers_predictions = esm_markers_predictions,
  rand_marker_predctions = rand_marker_predctions
)

save("all_marker_prediction_results", "marker_subset_prediction_results", 
     file = file.path(result_dir, "vp_stability_prediction.RData"))



