## S2MET Mapping
## Identify and map marker effect stability
## 
## This script will calculate marker effects across environment, test for significant
## marker effect stability, then attempt to map markers will stable or sensitive effects
## 

# Capture the trait to use (argument 1)
args <- commandArgs(trailingOnly = T)
tr <- args[1]

# List of packages to load
# List of packages
packages <- c("dplyr", "purrr", "tibble", "tidyr", "readr", "stringr", "readxl", "modelr", 
              "parallel", "purrrlyr", "rrBLUP", "ggplot2", "broom", "Matrix", "lme4qtl")

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

bopa_geno_dir <- geno_dir <-  "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/Data/Genos"
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

# Load an entry file
entry_list <- read_excel(file.path(entry_dir, "S2MET_project_entries.xlsx"))


# Grab the entry names that are not checks
tp <- entry_list %>% 
  filter(Class == "S2TP") %>% 
  pull(Line)

tp_geno <- intersect(tp, row.names(S2TP_imputed_multi_genos_mat))

# Define the checks
checks <- entry_list %>% 
  filter(Class == "Check") %>% 
  pull(Line)

entries <- entry_list %>% 
  pull(Line)

# Number of cores
n_cores <- detectCores()


# Matrix of genotype data for all individuals
# This must remain intact, and then will be subsetted below
# M <- s2_imputed_mat[c(tp_geno, vp_geno),]
M <- S2TP_imputed_multi_genos_mat[tp_geno,]


# Remove the environments in which the vp was only observed
S2_MET_BLUEs_use <- S2_MET_BLUEs %>%
  filter(line_name %in% c(tp_geno)) %>%
  mutate(line_name = as.factor(line_name),
         environment = as.factor(environment))


## Use the GWAS G model to estimate the effect of each marker in each environment
## This will be the environment-specific marker effect + the mean

# Subset markers by chromosome
markers_by_chrom <- S2TP_imputed_multi_genos_hmp %>% 
  select(marker = rs, chrom)

# All marker names
markers <- colnames(M)

# Create relationship matrices per chromosome
K_chr <- markers_by_chrom %>%
  split(.$chrom) %>%
  map(~M[,setdiff(markers, .$marker),drop = FALSE]) %>%
  map(~A.mat(X = ., min.MAF = 0, max.missing = 1))

## One trait at a time
trait_df <- S2_MET_BLUEs_use %>%
  filter(trait == tr) %>%
  droplevels()


  
# Model matrices for line_name and environment
Zg <- model.matrix(~ -1 + line_name, trait_df)
Ze <- model.matrix(~ -1 + environment, trait_df)

# Extract weights
wts <- trait_df$std_error^2

# Split markers by core
core_list <- markers_by_chrom %>% 
  mutate(core = sort(rep(seq(n_cores), length.out = nrow(.)))) %>%
  split(.$core)

# Iterate over the core list
marker_score_out <- mclapply(X = core_list, FUN = function(core) {
  
  # empty list to store results
  core_list_out <- vector("list", nrow(core)) %>%
    set_names(core$marker)
  
  # Apply a function over the marker matrix
  for (i in seq(nrow(core))) {
    # Subset the snp
    snp <- core[i,]
      
    mar <- Zg %*% M[,snp$marker, drop = FALSE]
    
    K_use <- K_chr[[snp$chrom]]
    
    # fit the model
    fit <- relmatLmer(value ~ mar:environment + (1|line_name) + (1|environment), trait_df, 
                      relmat = list(line_name = K_use), weights = wts)
    
    # Extract the coefficients
    core_list_out[[i]] <- tidy(fit) %>% 
      subset(str_detect(term, "mar"), c(term, estimate))
    
  }
  
  # return the list
  return(core_list_out)
  
}, mc.cores = n_cores)


# Save the output
save_file <- file.path(result_dir, str_c("S2MET_marker_by_env_", tr, ".RData"))
save("marker_score_out", file = save_file)



# ## Load the phenotype FW regression results
# load(file.path(result_dir, "S2MET_pheno_mean_fw_results.RData" ))
# 
# 
# ## Calculate marker effects in each environment for each trait
# S2_MET_marker_effect_env <- S2_MET_pheno_mean_fw %>% 
#   select(trait, line_name, environment, value, std_error) %>%
#   mutate(line_name = as.factor(line_name)) %>%
#   group_by(trait, environment) %>%
#   do({
#     # Extract the data.frame
#     df <- .
# 
#     # Model frame
#     mf <- model.frame(value ~ line_name, df)
#     # Vector of responses
#     y <- model.response(mf)
#     # Matrix of line names to be used to subset the genotypes
#     Zline <- model.matrix(~ -1 + line_name, mf)
# 
#     # Subset the genotypes
#     M1 <- Zline %*% M
#     
#     # Other matrices to use
#     X <- model.matrix(~ 1, mf)
#     K <- diag(ncol(M1))
#     
#     # Extract the weights for the R matrix
#     R <- diag(df$std_error^2)
#     
#     # Fit the model
#     # Results from this model correlate perfectly with the 'mixed.solve' model, so
#     # only using that model
#     # fit <- sommer::mmer(Y = y, X = X, Z = list(snp = list(Z = M1, K = K)), R = list(res = R))
# 
#     # Fit the mixed model
#     fit <- mixed.solve(y = y, Z = M1)
# 
#     # Extract the marker effects and convert to data.frame, then output
#     fit$u %>%
#       data_frame(marker = names(.), effect = ., beta = fit$beta) })
# 
# # Ungroup
# S2_MET_marker_effect_env <- ungroup(S2_MET_marker_effect_env)
# 
# ### Perform FW regression using the marker effects
# 
# 
# 
# ## Combine the marker effects by environment with the environmental effect from
# ## FW regression
# S2_MET_marker_env_eff <- left_join(distinct(S2_MET_pheno_mean_fw, trait, environment, h),
#                                    S2_MET_marker_effect_env, by = c("environment", "trait"))
# 
# 
# # Map over the list and calculate FW regression slopes and deviations for each
# # of the markers
# marker_stability_coef <- S2_MET_marker_env_eff %>%
#   group_by(trait, marker) %>%
#   do({
#     # Extract the data
#     df <- .
#     
#     # # Add the beta term to the effect
#     # df1 <- df %>%
#     #   mutate(mar_value = effect + beta)
# 
#     # Run the regression and tidy the results
#     fit <- lm(effect ~ h, data = df)
#     # fit <- lm(mar_value ~ h, data = df1)
#     fit_tidy <- tidy(fit)
# 
#     # Return coefficients
#     data.frame(
#       b = coef(fit)[2],
#       b_std_error = subset(fit_tidy, term == "h", std.error, drop = T),
#       df = df.residual(fit),
#       delta = mean(resid(fit)^2),
#       row.names = NULL
#     ) })
# 
# 
# ## Add marker information (position, etc.)
# marker_stability_coef_info <- marker_stability_coef %>%
#   left_join(select(S2TP_imputed_multi_genos_hmp, marker = rs, alleles:cM_pos), ., by = "marker")
# 
# # Add the regression coefficients to the marker effects
# S2_MET_marker_eff_pheno_fw <- left_join(S2_MET_marker_env_eff, marker_stability_coef_info,
#                                         by = c("marker", "trait")) %>%
#   # Rearrange some columns and edit some variables
#   select(., trait, environment, env_effect = h, marker, chrom, alleles, pos, cM_pos,
#          mar_effect = effect, b_std_error, df, stability_term, estimate)
# 
# 
# ## Save the results
# save_file <- file.path(result_dir, "S2MET_pheno_mar_eff_fw_results.RData")
# save("S2_MET_marker_effect_env", "S2_MET_marker_eff_pheno_fw", file = save_file)

# n_cores <- detectCores()
# n_perm <- 100

# # Perform permutation analysis by shuffling the values of traits in each environment
# phenos_mxe_perm <- S2_MET_BLUEs_use %>%
#   mutate(environment = as.character(environment)) %>%
#   split(.$trait) %>% 
#   map(~split(., .$environment) %>% 
#         map(~permute(data = ., n = n_perm, line_name)) %>% 
#         map("perm") %>%
#         transpose() %>%
#         map(., ~map_df(., as_data_frame)))
# 
# 
# # Map over the traits and permutation to calculate environment-specific marker
# # effects
# 
# # Split the iterations for parallelization
# phenos_mxe_perm_split <- phenos_mxe_perm %>% 
#   transpose() %>% 
#   map(bind_rows) %>%
#   split(., cut(seq_along(.), breaks = n_cores))
# 
# phenos_mxe_perm_parallel_out <- mclapply(X = phenos_mxe_perm_split, function(core) {
#   
#   map(core, function(iter) {
#     
#     out <- group_by(iter, trait, environment) %>%
#       do({
#         # Extract the data
#         df <- .
#         
#         mf <- model.frame(value ~ line_name, df)
#         y <- model.response(mf)
#         Zline <- model.matrix(~ -1 + line_name, mf)
#         # Subset the genotypes
#         M1 <- Zline %*% M
#         # Fit the mixed model
#         fit <- mixed.solve(y = y, Z = M1)
#         
#         fit$u %>% 
#           data.frame(marker = names(.), effect = ., beta = fit$beta, 
#                      row.names = NULL, stringsAsFactors = FALSE)
#         
#       }) })
#   
#   map(out, ungroup) }, mc.cores = n_cores)
#       
# 
# # Save the data
# save_file <- file.path(result_dir, "S2MET_pheno_mar_eff_by_env_perm.RData")
# save("phenos_mxe_perm_parallel_out", file = save_file)


 
# ## Calculate marker effects x environment as in Lopez-Cruz 2015
# 
# ## Group by trait and estimate the base marker effect and enivornmnet-specific deviation
# mar_effect_by_env <- S2_MET_BLUEs_use %>%
#   mutate(line_name = as.factor(line_name),
#          environment = as.factor(environment)) %>%
#   group_by(trait) %>%
#   do({
# 
#     # Convert the . to df
#     df <- .
#     # Create a model frame with lines and environments
#     mf <- model.frame(value ~ line_name + environment, data = df, drop.unused.levels = TRUE)
#     env_names <- levels(mf$environment)
# 
#     # Extract the marker names
#     mar_names <- colnames(M)
#     
#     # Response vector
#     y <- model.response(mf)
# 
#     # Model matrix for the fixed effects
#     X <- model.matrix(~ environment, data = mf)
# 
#     # Model matrix of the line_name incidence
#     Z_line <- sparse.model.matrix(~ -1 + line_name, mf)
#     # Model matrix of the marker random effects
# 
#     ## First the main effect of each marker
#     Z0 <- Z_line %*% M
#     K0 <- Diagonal(ncol(Z0))
# 
#     # Random deviation of each marker in each environment
#     Z1_list <- mf %>%
#       split(.$environment) %>%
#       map(~ sparse.model.matrix(~ -1 + line_name, .)) %>%
#       map(~ . %*% M)
#     
#     Z1 <- .bdiag(Z1_list)
#     colnames(Z1) <- str_c(rep(mar_names, length(env_names)), 
#                           rep(env_names, each = length(mar_names)), sep = ":")
#     
#     K1 <- Diagonal(ncol(Z1))
# 
#     # fit <- EMMREML::emmremlMultiKernel(y = y, X = X, Zlist = list(Z0, Z1), Klist = list(K0, K1))
#     fit <- sommer::mmer(Y = y, X = X, Z = list(m = list(Z = Z0, K = K0), mxe = list(Z = Z1, K = K1)))
# 
#     m <- fit$u.hat$m %>% 
#       data.frame(marker = row.names(.), environment = NA, effect = ., 
#                  row.names = NULL, stringsAsFactors = FALSE) %>% 
#       dplyr::rename(effect = T1)
#     
#     # Rename the random effect vector
#     mxe <- fit$u.hat$mxe %>% 
#       data.frame(marker = row.names(.), effect = ., row.names = NULL, stringsAsFactors = FALSE) %>% 
#       separate(marker, c("marker", "environment"), sep = ":", fill = "right") %>% 
#       dplyr::rename(effect = T1)
#     
#     # Return
#     bind_rows(m, mxe)
# 
#   })
# 
# ## Save this
# save_file <- file.path(result_dir, "S2MET_marker_eff_by_env.RData")
# save("mar_effect_by_env", file = save_file)




# # Map over the list and calculate FW regression slopes and deviations for each
# # of the markers
# marker_stability_coef <- S2_MET_marker_env_eff %>%
#   map(~group_by(., trait, marker) %>%
#         do({
#           # Extract the data
#           df <- .
#           
#           # Run the regression and tidy the results
#           fit <- lm(effect ~ h, data = df)
#           fit_tidy <- tidy(fit)
#           
#           # Return coefficients
#           data.frame(
#             b = subset(fit_tidy, term == "h", estimate, drop = T),
#             b_std_error = subset(fit_tidy, term == "h", std.error, drop = T),
#             df = df.residual(fit),
#             delta = mean(resid(fit)^2),
#             row.names = NULL
#           ) }) %>%
#         gather(stability_term, estimate, -trait, -marker, -b_std_error, -df) %>%
#         ungroup() )
# 
# 
# ## Add marker information (position, etc.)
# marker_stability_coef_info <- marker_stability_coef %>% 
#   map(~left_join(rename(snp_info, marker = `rs#`), ., by = "marker"))
# 
# # Add the regression coefficients to the marker effects
# S2_MET_marker_eff_pheno_fw <- list(S2_MET_marker_env_eff, marker_stability_coef_info) %>%
#   pmap(~left_join(.x, .y, by = c("marker", "trait"))) %>%
#   # Rearrange some columns and edit some variables
#   map(~select(., trait, environment, env_effect = h, marker, chrom, alleles, pos, cM_pos, 
#               mar_effect = effect, b_std_error, df, stability_term, estimate))


