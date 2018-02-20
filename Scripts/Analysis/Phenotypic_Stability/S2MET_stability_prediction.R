## Genomic prediction of FW sensitivity coefficients
## 
## Conduct genomic prediction and cross-validation of FW sensitivity coefficients
## For CV, use both GBS markers and BOPA markers. Also conduct predictions using the 
## genotypic means
## 

# List of packages to load
# List of packages
packages <- c("dplyr", "purrr", "tibble", "tidyr", "readr", "stringr", "readxl", "modelr", 
              "parallel", "purrrlyr", "rrBLUP")

# Set the directory of the R packages
package_dir <- NULL
package_dir <- "/panfs/roc/groups/6/smithkp/neyha001/R/x86_64-pc-linux-gnu-library/3.4/"

# Load all packages
invisible(lapply(packages, library, character.only = TRUE, lib.loc = package_dir))


## Directories
proj_dir <- "C:/Users/Jeff/Google Drive/Barley Lab/Projects/S2MET_Mapping//"
proj_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET_Mapping/" 

alt_proj_dir <- "C:/Users/Jeff/Google Drive/Barley Lab/Projects/S2MET"
alt_proj_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET/"

# Geno, pheno, and enviro data
geno_dir <-  "C:/Users/Jeff/Google Drive/Barley Lab/Projects/Genomics/Genotypic_Data/GBS_Genotype_Data/"
bopa_geno_dir <- "C:/Users/Jeff/Google Drive/Barley Lab/Projects/Genomics/Genotypic_Data/BOPA_Genotype_Data/"
pheno_dir <- file.path(alt_proj_dir, "Phenotype_Data/")
env_var_dir <- file.path(alt_proj_dir, "Environmental_Variables")

geno_dir <- bopa_geno_dir <-  "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/Data/Genos"
pheno_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/Data/Phenos"
env_var_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/Data/Environmental_Data"


# Other directories
fig_dir <- file.path(proj_dir, "Figures/")
map_dir <- file.path(proj_dir, "Mapping")
entry_dir <- file.path(alt_proj_dir, "Plant_Materials")
analysis_dir <- file.path(proj_dir, "Analysis")
result_dir <- file.path(proj_dir, "Results")

# Load the geno data
load(file.path(bopa_geno_dir, "S2TP_multi_genos.RData"))

# Load the FW results
load(file.path(result_dir, "S2MET_pheno_mean_fw_results.RData"))
# Load the marker FW results
load(file.path(result_dir, "S2MET_marker_mean_fw_results.RData"))

# Find the most stable/sensitive and average markers
S2_MET_marker_mean_fw_tidy <- S2_MET_marker_mean_fw %>% 
  distinct(marker, chrom, pos, trait, b, delta) %>%
  mutate(log_delta = log(delta), b = b + 1) %>%
  select(-delta) %>% 
  gather(coef, estimate, -marker:-trait)

# What should be the significance level
alpha <- 0.05

# For each trait, calculate empirical thresholds for significance
S2_MET_marker_fw_sig <- S2_MET_marker_mean_fw_tidy %>%
  filter(coef == "b") %>% 
  group_by(trait) %>% 
  # mutate(estimate = scale(estimate)) %>%
  mutate(lower_perc = quantile(estimate, alpha / 2), 
         upper_perc = quantile(estimate, 1 - (alpha / 2))) %>%
  ungroup() %>%
  mutate(significance = case_when(estimate >= upper_perc ~ "sensitive",
                                  estimate <= lower_perc ~ "stable",
                                  TRUE ~ "average"),
         marker_type = if_else(str_detect(marker, "^S"), "GBS", "BOPA"))

# Load an entry file
entry_list <- read_excel(file.path(entry_dir, "S2MET_project_entries.xlsx"))


# Grab the entry names that are not checks
tp <- entry_list %>% 
  filter(Class == "S2TP") %>% 
  pull(Line)

# Find the tp and vp that are genotypes
tp_geno <- intersect(tp, row.names(S2TP_imputed_multi_genos_mat))

# Rename the marker matrix
M <- S2TP_imputed_multi_genos_mat[tp_geno,]

# Overall K
K <- A.mat(X = M, min.MAF = 0, max.missing = 1)


# Create a tidy dataset to model
S2_MET_pheno_mean_fw_tomodel <- S2_MET_pheno_mean_fw %>% 
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
cv_parts <- S2_MET_pheno_mean_fw_tomodel %>%
  group_by(trait, coef) %>%
  do(folds = rerun(.n = n_cv_iter, crossv_kfold(data = ., k = n_k)))

# Run the cross-validation
cv_results_all <- cv_parts %>%
  do(map_df(.$folds, ~rowwise(.) %>% do(predict_RR(train = .$train, test = .$test, K = K)) %>% 
              ungroup() %>% mutate(acc = cor(value, pred_value)) %>%
              distinct(trait, coef,  acc)))
              


# Now use different marker sets
marker_stability_sets <- S2_MET_marker_fw_sig %>% 
  group_by(trait, significance) %>% 
  do(K = A.mat(X = M[,.$marker, drop = FALSE], min.MAF = 0, max.missing = 1))

# Combine
cv_parts_marker_sets <- full_join(cv_parts, marker_stability_sets, by = "trait")

# Run the cross-validation
cv_results_marker_sets <- cv_parts_marker_sets %>%
  group_by(trait, coef, significance) %>%
  do({
    flds <- .$folds[[1]]
    K_use <- .$K[[1]]
    
    map_df(flds, ~rowwise(.) %>% do(predict_RR(train = .$train, test = .$test, K = K_use)) %>% 
              ungroup() %>% mutate(acc = cor(value, pred_value)) %>%
              distinct(trait, coef,  acc)) })


# Combine the data.frames
cv_results <- bind_rows(
  ungroup(cv_results_all) %>% mutate(significance = "all_markers"),
  ungroup(cv_results_marker_sets)
)



## Randomly subset the same number of markers for creating relationship matrices
# Then run 50 k-fold CV replications per marker sample - average
# Repeat this 100 times
n_cv_iter <- 10
n_samples <- 100


# First nest the marker groups
S2_MET_marker_fw_nest <- S2_MET_marker_fw_sig %>% 
  group_by(trait, significance) %>% 
  nest(marker) %>%
  bind_rows(., S2_MET_marker_fw_sig %>% group_by(trait) %>% 
              nest(marker) %>% mutate(significance = "all"))

# Find the minimum number of markers
min_markers <- S2_MET_marker_fw_nest$data %>% 
  map_dbl(~nrow(.)) %>% 
  min()

# Make random samples of the markers
S2_MET_marker_fw_samples <- S2_MET_marker_fw_nest %>% 
  group_by(trait, significance) %>%
  do(data_frame(rep = seq(n_samples), sample = rerun(.n = n_samples, sample_n(.$data[[1]], size = min_markers))))

# Create relationship matrices
S2_MET_marker_fw_samples_K <- S2_MET_marker_fw_samples %>% 
  ungroup() %>%
  mutate(K = map(sample, ~A.mat(X = M[,.$marker, drop = FALSE], min.MAF = 0, max.missing = 1)))

# Combine with the phenotypic data
S2_MET_marker_fw_samples_tomodel <- left_join(
  S2_MET_marker_fw_samples_K,
  S2_MET_pheno_mean_fw_tomodel %>% group_by(trait, coef) %>% nest(),
  by = "trait"
)

## Create cv reps for each line and run the CV, report the mean
# Make an empty list
cv_results_samples <- vector("list", nrow(S2_MET_marker_fw_samples_tomodel))

# Iterate over the rows
for (i in seq(nrow(S2_MET_marker_fw_samples_tomodel))) {
  
  df <- S2_MET_marker_fw_samples_tomodel[i,]
  
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
save_file <- file.path(result_dir, "S2MET_stability_crossv_results.RData")
save("cv_results", "cv_results_samples_tidy", file = save_file)

