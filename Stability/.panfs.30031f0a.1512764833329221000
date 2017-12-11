## Find the proportion of GxE variance explained by stable and non-significant markers
## 
## 

# List of packages to load
# List of packages
packages <- c("dplyr", "purrr", "tibble", "tidyr", "readr", "stringr", "readxl", "parallel", "rrBLUP",
              "purrrlyr")

# Set the directory of the R packages
package_dir <- NULL
package_dir <- "/panfs/roc/groups/6/smithkp/neyha001/R/x86_64-pc-linux-gnu-library/3.4/"
# package_dir <- "/panfs/roc/groups/6/smithkp/neyha001/R/x86_64-unknown-linux-gnu-library/3.2/"

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


# Load the phenotypic data
load(file.path(pheno_dir, "S2_MET_BLUEs.RData"))
# Load the genotypic data
# load(file.path(geno_dir, "S2_genos_mat.RData"))
# load(file.path(geno_dir, "S2_genos_hmp.RData"))
load(file.path(bopa_geno_dir, "S2TP_multi_genos.RData"))
# Load environmental data
load(file.path(env_var_dir, "environmental_data_compiled.RData"))

# Load the marker effect stability results
load(file.path(result_dir, "S2MET_pheno_mar_eff_fw_results.RData"))

# Load an entry file
entry_list <- read_excel(file.path(entry_dir, "S2MET_project_entries.xlsx"))


# Grab the TP line names
tp <- entry_list %>% 
  filter(Class == "S2TP") %>% 
  pull(Line)

# Find the intersection of TP lines and those that have been genotyped
tp_geno <- intersect(tp, row.names(S2TP_imputed_multi_genos_mat))


# Filter the BLUEs to use
S2_MET_BLUEs_use <- S2_MET_BLUEs %>% 
  filter(line_name %in% c(tp_geno)) %>%
  droplevels()

# Create a genotype marker matrix
M <- S2TP_imputed_multi_genos_mat[tp_geno,]


## First group markers into those that are highly stable, highly sensitive, or neither


# Get the distinct observations of the stability measures
S2_MET_marker_eff_pheno_fw_uniq <- S2_MET_marker_eff_pheno_fw %>% 
  select(trait, marker:pos, b_std_error, df, stability_term, estimate) %>% 
  distinct()

# What should be the significance level
alpha <- 0.05

# For each trait, calculate empirical thresholds for significance
S2_MET_marker_eff_pheno_fw_sig <- S2_MET_marker_eff_pheno_fw_uniq %>%
  filter(stability_term == "b") %>% 
  group_by(trait) %>% 
  # mutate(estimate = scale(estimate)) %>%
  mutate(lower_perc = quantile(estimate, alpha / 2), 
         upper_perc = quantile(estimate, 1 - (alpha / 2))) %>%
  ungroup() %>%
  mutate(significance = case_when(estimate >= upper_perc ~ "sensitive",
                                  estimate <= lower_perc ~ "stable",
                                  TRUE ~ "not_significant"),
         marker_type = if_else(str_detect(marker, "^S"), "GBS", "BOPA"))


## Fit the mixed models to determine variance components

# Number of model fittings
n_iter <- 10
# Detect cores
n_cores <- detectCores()



## Calculate relationship matrices
# Main effect relationship matrix
K <- A.mat(X = M, min.MAF = 0, max.missing = 1)

# Extract the traits
trts <- unique(S2_MET_marker_eff_pheno_fw_sig$trait)
# Create a results list
var_comp_list <- vector("list", length(trts)) %>%
  set_names(trts)

# Iterate over the traits
for (tr in trts) {
  
  # Extract the data to use
  pheno_df <- S2_MET_BLUEs_use %>%
    filter(trait == tr) %>%
    # filter(environment %in% sample(unique(.$environment), 5)) %>%
    droplevels()
  
  # Extract the significant markers
  sig_mar_df <- S2_MET_marker_eff_pheno_fw_sig %>% 
    filter(trait == tr, significance != "not_significant")
  # Extract non-sig markers
  non_sig_mar_df <- S2_MET_marker_eff_pheno_fw_sig %>% 
    filter(trait == tr, significance == "not_significant")
  
  ## Subset marker relationship matrices
  sig_mar <- pull(sig_mar_df, marker)
  M_s <- M[,sig_mar]
  K_s <- A.mat(X = M_s, min.MAF = 0, max.missing = 1)

  non_sig_mar <- pull(non_sig_mar_df, marker)
  # Sample markers to construct Ks
  non_sig_mar_sample <- replicate(n = n_iter, {sort(sample(non_sig_mar, size = length(sig_mar)))}, simplify = FALSE)
  M_ns_list <- map(non_sig_mar_sample, ~M[,.])
  K_ns_list <- map(M_ns_list, ~A.mat(X = ., min.MAF = 0, max.missing = 1))
  
  ## Model frame/matrices
  mf <- model.frame(value ~ line_name + environment + std_error, pheno_df) %>%
    mutate(gt = interaction(line_name, environment)) %>%
    droplevels()
  
  y <- mf$value
  # Fixed effect of environments?
  X <- model.matrix(~ environment, mf)
  
  # Random effect of environment
  Z_e <- model.matrix(~ -1 + environment, mf)
  
  # Random effect of genotypic main effect
  Z_g <- model.matrix(~ -1 + line_name, mf)
  
  # Random effect of GxE
  Z_gt <- model.matrix(~ -1 + gt, mf, drop.unused.levels = T)
  
  ## K_s used in the model
  K_s_use <- sommer::hadamard.prod(x = (Z_g %*% K_s %*% t(Z_g)), y = tcrossprod(Z_e))
  
  ## Split the list of K_ns by core
  K_ns_list_core <- split(x = K_ns_list, f = cut(x = seq(n_iter), breaks = n_cores)) %>%
    set_names(seq(n_cores))
  # K_ns_list_core <- list(`1` = K_ns_list)
  
  # Parallelize over the core split list of K_ns matrices
  parallel_out <- mclapply(X = K_ns_list_core, FUN = function(K_ns_core) {
    
    # Calculate the K matrices that will be used
    K_ns_core_use <- vector("list", length(K_ns_core))
    
    for (i in seq_along(K_ns_core_use)) {
      K_ns_core_use[[i]] <- sommer::hadamard.prod(x = (Z_g %*% K_ns_core[[i]] %*% t(Z_g)), y = tcrossprod(Z_e))
    }
    
    
    # Map over the list
    var_comp_out <- K_ns_core_use %>%
      map_df(~sommer::mmer(Y = y, X = X, Z = list(g = list(Z = Z_g, K = K), gt_s = list(Z = Z_gt, K = K_s_use),
                                               gt_ns = list(Z = Z_gt, K = .))) %>% .$var.comp %>% map_df(unname))
    
    return(var_comp_out)
    
  }, mc.cores = n_cores)
  
  # Add the output to the list
  var_comp_list[[tr]] <- parallel_out
  
} # Close the trait for loop

# Save the data
save_file <- file.path(result_dir, "S2MET_pheno_mar_fw_varcomp.RData")
save("var_comp_list", file = save_file)
