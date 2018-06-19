## Genomic prediction of FW sensitivity coefficients
## 
## Conduct genomic prediction and cross-validation of FW sensitivity coefficients
## For CV, use both GBS markers and BOPA markers. Also conduct predictions using the 
## genotypic means
## 
## Author: Jeff Neyhart
## Last updated: June 18, 2018
## 


# Run the source script - MSI
repo_dir <- "/panfs/roc/groups/6/smithkp/neyha001/QTLMapping/S2MET_Mapping/"
source(file.path(repo_dir, "source_MSI.R"))

# # Run the source script - local
# repo_dir <- getwd()
# source(file.path(repo_dir, "source.R"))

# Load the FW results
load(file.path(result_dir, "pheno_mean_fw_results.RData"))
# Load the marker FW results
load(file.path(result_dir, "marker_mean_fw_results.RData"))


# Find the most plastic markers and stable markers
marker_mean_fw_tidy <- marker_mean_fw %>% 
  distinct(marker, chrom, pos, trait, b, delta) %>%
  mutate(log_delta = log(delta), b = b + 1) %>%
  select(-delta) %>% 
  gather(coef, estimate, -marker:-trait)

# What should be the significance level
alpha <- 0.05

# For each trait, calculate empirical thresholds for significance
marker_fw_sig <- marker_mean_fw_tidy %>%
  filter(coef == "b") %>% 
  group_by(trait) %>% 
  # mutate(estimate = scale(estimate)) %>%
  mutate(lower_perc = quantile(estimate, alpha / 2), 
         upper_perc = quantile(estimate, 1 - (alpha / 2))) %>%
  ungroup() %>%
  mutate(significance = case_when(estimate >= upper_perc ~ "plastic",
                                  estimate <= lower_perc ~ "plastic",
                                  TRUE ~ "stable"),
         marker_type = if_else(str_detect(marker, "^S"), "GBS", "BOPA"))



# Rename the marker matrix
M <- S2TP_imputed_multi_genos_mat[tp_geno,]
# Overall K
K <- A.mat(X = M, min.MAF = 0, max.missing = 1)


# Create a tidy dataset to model
pheno_mean_fw_tomodel <- S2_MET_pheno_mean_fw %>% 
  distinct(trait, line_name, g, b, delta) %>%
  mutate(line_name = as.factor(line_name),
         log_delta = log(delta)) %>%
  select(-delta) %>%
  gather(coef, value, g:log_delta)

# Write a function that takes a train and test set and predicts using rrBLUP
predict_RR <- function(train, test, K) {
  
  # Convert to df
  train_df <- as.data.frame(train)
  test_df <- as.data.frame(test)
  
  # Create the model matrix
  mf <- model.frame(value ~ line_name, train_df)
  y <- model.response(mf)
  Z <- model.matrix(~ -1 + line_name, mf)
  
  fit <- mixed.solve(y = y, Z = Z, K = K)
  
  # Tidy
  u_hat_tidy <- fit$u %>% 
    data.frame(line_name = names(.), pred_value = ., stringsAsFactors = FALSE, row.names = NULL)
  
  # Combine and return the predictions
  suppressWarnings(left_join(test_df, u_hat_tidy, by = "line_name"))
  
}


# Set the number of CV iterations and the number of K folds
n_cv_iter <- 100
n_k <- 7

# Create the cv parts
cv_parts <- pheno_mean_fw_tomodel %>%
  group_by(trait, coef) %>%
  do(folds = rerun(.n = n_cv_iter, crossv_kfold(data = ., k = n_k)))

# Run the cross-validation
cv_results_all <- cv_parts %>%
  do(map_df(.$folds, ~rowwise(.) %>% do(predict_RR(train = .$train, test = .$test, K = K)) %>% 
              ungroup() %>% mutate(acc = cor(value, pred_value)) %>%
              distinct(trait, coef,  acc)))
              


## Randomly subset the same number of markers for creating relationship matrices
# Then run 10 k-fold CV replications per marker sample - average
# Repeat this 100 times
n_cv_iter <- 10
n_samples <- 100


# First nest the marker groups
marker_fw_nest <- marker_fw_sig %>% 
  group_by(trait, significance) %>% 
  nest(marker) %>%
  bind_rows(., marker_fw_sig %>% group_by(trait) %>% 
              nest(marker) %>% mutate(significance = "all"))

# Combine the stable and sensitive markers into a group called "mix"
marker_mix_nest <- marker_fw_nest %>% 
  filter(significance %in% c("stable", "sensitive")) %>% 
  group_by(trait) %>% 
  do(data = bind_rows(.$data)) %>%
  ungroup() %>%
  mutate(significance = "mixed")

# Combine
marker_fw_nest_use <- bind_rows(marker_fw_nest, marker_mix_nest)

# Find the minimum number of markers
min_markers <- marker_fw_nest_use$data %>% 
  map_dbl(~nrow(.)) %>% 
  min()

# Make random samples of the markers
marker_fw_samples <- marker_fw_nest_use %>% 
  group_by(trait, significance) %>%
  do(data_frame(rep = seq(n_samples), sample = rerun(.n = n_samples, sample_n(.$data[[1]], size = min_markers))))

# Create relationship matrices
marker_fw_samples_K <- marker_fw_samples %>% 
  ungroup() %>%
  mutate(K = map(sample, ~A.mat(X = M[,.$marker, drop = FALSE], min.MAF = 0, max.missing = 1)))

# Combine with the phenotypic data
marker_fw_samples_tomodel <- left_join(
  marker_fw_samples_K,
  pheno_mean_fw_tomodel %>% group_by(trait, coef) %>% nest(),
  by = "trait"
)

## Create cv reps for each line and run the CV, report the mean
# Make an empty list
cv_results_samples <- vector("list", nrow(marker_fw_samples_tomodel))

# Iterate over the rows
for (i in seq(nrow(marker_fw_samples_tomodel))) {
  
  df <- marker_fw_samples_tomodel[i,]
  
  trait <- df$trait
  coef <- df$coef
  rep <- df$rep
  significance <- df$significance
  
  # Create CV reps
  cv_parts <- data_frame(
    iter = seq(n_cv_iter),
    parts = rerun(.n = n_cv_iter, crossv_kfold(data = df$data[[1]], k = n_k))
  ) %>% unnest()
  
  # Map over them and calculate cv accuracy
  cv_out <- cv_parts %>% 
    group_by(iter, .id) %>%
    do(cv_out = predict_RR(train = .$train[[1]], test = .$test[[1]], K = df$K[[1]]))
  
  # Return the results
  acc <- cv_out %>% 
    unnest(cv_out) %>% 
    group_by(iter) %>% 
    summarize(rep_acc = cor(value, pred_value))
  
  # Add to the list
  cv_results_samples[[i]] <- list(
    trait = trait,
    coef = coef,
    rep = rep,
    significance = significance,
    acc = acc
  )
  
}

# Tidy up
cv_results_samples_tidy <- cv_results_samples %>%
  map_df(as.data.frame) %>%
  tbl_df()

# Save this
save_file <- file.path(result_dir, "stability_crossv_results.RData")
save("cv_results", "cv_results_samples_tidy", file = save_file)

