## S2MET Mapping
## Identify and map marker effect stability
## 
## This script will calculate marker effects across environment, test for significant
## marker effect stability, then attempt to map markers will stable or sensitive effects
## 

# List of packages to load
# List of packages
packages <- c("dplyr", "purrr", "tibble", "tidyr", "readr", "stringr", "readxl", 
              "parallel", "purrrlyr", "rrBLUP", "FW", "ggplot2", "broom", "Matrix")

# Set the directory of the R packages
package_dir <- NULL
# package_dir <- "/panfs/roc/groups/6/smithkp/neyha001/R/x86_64-pc-linux-gnu-library/3.4/"

# Load all packages
invisible(lapply(packages, library, character.only = TRUE, lib.loc = package_dir))


## Directories
proj_dir <- "C:/Users/Jeff/Google Drive/Barley Lab/Projects/S2MET_Mapping//"
# proj_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET_Mapping/" 

alt_proj_dir <- "C:/Users/Jeff/Google Drive/Barley Lab/Projects/S2MET"
# alt_proj_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET/"

# Geno, pheno, and enviro data
geno_dir <-  "C:/Users/Jeff/Google Drive/Barley Lab/Projects/Genomics/Genotypic_Data/GBS_Genotype_Data/"
pheno_dir <- file.path(alt_proj_dir, "Phenotype_Data/")
env_var_dir <- file.path(alt_proj_dir, "Environmental_Variables")

# geno_dir <-  "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/Data/GBS_Genos"
# pheno_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/Data/Phenos"
# env_var_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/Data/Environmental_Data"


# Other directories
fig_dir <- file.path(proj_dir, "Figures/")
map_dir <- file.path(proj_dir, "Mapping")
entry_dir <- file.path(alt_proj_dir, "Plant_Materials")
analysis_dir <- file.path(proj_dir, "Analysis")
result_dir <- file.path(proj_dir, "Results")


# Load the phenotypic data
load(file.path(pheno_dir, "S2_MET_BLUEs.RData"))
# Load the genotypic data
load(file.path(geno_dir, "S2_genos_mat.RData"))
# load(file.path(geno_dir, "S2_genos_hmp.RData"))
# Load environmental data
load(file.path(env_var_dir, "environmental_data_compiled.RData"))

# Load an entry file
entry_list <- read_excel(file.path(entry_dir, "S2MET_project_entries.xlsx"))


# Grab the entry names that are not checks
tp <- entry_list %>% 
  filter(Class == "S2TP") %>% 
  pull(Line)

vp <- entry_list %>% 
  filter(Class == "S2C1R") %>% 
  pull(Line)

# Find the tp and vp that are genotypes
tp_geno <- intersect(tp, row.names(s2_imputed_mat))
vp_geno <- intersect(vp, row.names(s2_imputed_mat))

# Define the checks
checks <- entry_list %>% 
  filter(Class == "Check") %>% 
  pull(Line)

entries <- entry_list %>% 
  pull(Line)

# Define a list of two populations
pop_list <- list(pop_all = c(tp_geno, vp_geno),
                 pop_tp = tp_geno)

# Matrix of genotype data for all individuals
# This must remain intact, and then will be subsetted below
M <- s2_imputed_mat[c(tp_geno, vp_geno),]


# Remove the environments in which the vp was only observed
S2_MET_BLUEs_use <- S2_MET_BLUEs %>%
  # filter(!grepl(pattern = "BZI|HTM", x = environment)) %>% # Let this in.
  group_by(trait, environment) %>%
  filter(sum(line_name %in% tp) > 1,
         line_name %in% c(tp_geno, vp_geno)) %>%
  ungroup() %>%
  mutate(line_name = as.factor(line_name),
         environment = as.factor(environment))


## Calculate marker effects in each environment for each trait
## Do this for the 'all' population and the 'tp' population.
S2_MET_marker_effect_env_list <- pop_list %>%
  map(function(pop) {
    S2_MET_BLUEs_use %>% 
      filter(line_name %in% pop) %>% 
      group_by(trait, environment) %>% 
      do({
        # Extract the data.frame
        df <- .
        
        # Model frame
        mf <- model.frame(value ~ line_name, df)
        # Vector of responses
        y <- model.response(mf)
        # Matrix of line names to be used to subset the genotypes
        Zline <- model.matrix(~ -1 + line_name, mf)
        
        # Subset the genotypes
        M1 <- Zline %*% M
          
        # Fit the mixed model
        fit <- mixed.solve(y = y, Z = M1)
        
        # Extract the marker effects and convert to data.frame, then output
        fit$u %>% 
          data_frame(marker = names(.), effect = ., beta = fit$beta) }) %>%
      ungroup() })

# Save
save_file <- file.path(result_dir, "S2MET_marker_eff_by_env")
save("S2_MET_marker_effect_env_list", file = save_file)

# 
# 
# ## Calculate marker effects x environment as in Lopez-Cruz 2015
# # 
# # # Test on three environments
# # test_blues <- S2_MET_BLUEs_use %>%
# #   filter(trait == "GrainYield", environment %in% c("BCW16", "STP16", "CRM16"))
# # 
# # # Training data
# # train_blues <- test_blues %>%
# #   filter(line_name %in% tp_geno) %>%
# #   droplevels()
# 
# # Marker matrix
# M <- s2_imputed_mat_use
# 
# mar_names <- colnames(M)
# 
# 
# ## Group by trait and estimate the base marker effect and enivornmnet-specific deviation
# mar_effect_by_env <- S2_MET_BLUEs_use %>%
#   group_by(trait) %>%
#   do({
#     
#     # Convert the . to df
#     df <- .
#     # Create a model frame with lines and environments
#     mf <- model.frame(value ~ line_name + environment, data = df, drop.unused.levels = TRUE)
#     env_names <- levels(mf$environment)
# 
#     # Response vector
#     y <- model.response(mf)
# 
#     # Model matrix for the fixed effects
#     X <- model.matrix(~ -1 + environment, data = mf)
# 
#     # Model matrix of the line_name incidence
#     Z_line <- model.matrix(~ -1 + line_name, mf)
#     # Model matrix of the marker random effects
#     
#     ## First the main effect of each marker
#     Z0 <- Z_line %*% M
#     K0 <- Diagonal(ncol(Z0))
#     
#     # Random deviation of each marker in each environment
#     Z1 <- mf %>% 
#       split(.$environment) %>% 
#       map(~ model.matrix(~ -1 + line_name, .)) %>% 
#       map(~ . %*% M)
#     Z1 <- .bdiag(Z1)
#     
#     K1 <- Diagonal(ncol(Z1))
#     
#     fit <- EMMREML::emmremlMultiKernel(y = y, X = X, Zlist = list(Z0, Z1), Klist = list(K0, K1))
# 
#     # Rename the random effect vector
#     fit$uhat %>%
#       setNames(c(mar_names, paste(rep(mar_names, length(env_names)), 
#                                   rep(env_names, each = length(mar_names)), sep = ":"))) %>%
#       data.frame(marker = names(.), effect = ., row.names = NULL, stringsAsFactors = FALSE) %>%
#       separate(marker, c("marker", "environment"), sep = ":", fill = "right")
#     
#   })
# 
# ## Save this
# save_file <- file.path(result_dir, "S2MET_marker_eff_by_env.RData")
# save("mar_effect_by_env", file = save_file)



### Perform FW regression using the marker effects

## Load the phenotype FW regression results
load(file.path(result_dir, "S2MET_pheno_fw_regression_results.RData" ))

## Combine the marker effects by environment with the environmental effect from
## FW regression
S2_MET_marker_env_eff <- list(S2_MET_pheno_fw, S2_MET_marker_effect_env_list) %>% 
  pmap(~left_join(x = distinct(.x, trait, environment, h), y = .y, by = c("environment", "trait")))


# Map over the list and calculate FW regression slopes and deviations for each
# of the markers
marker_stability_coef <- S2_MET_marker_env_eff %>%
  map(~group_by(., trait, marker) %>%
        do({
          # Extract the data
          df <- .
          
          # Run the regression and tidy the results
          fit <- lm(effect ~ h, data = df)
          fit_tidy <- tidy(fit)
          
          # Return coefficients
          data.frame(
            b = subset(fit_tidy, term == "h", estimate, drop = T),
            b_std_error = subset(fit_tidy, term == "h", std.error, drop = T),
            df = df.residual(fit),
            delta = mean(resid(fit)^2),
            row.names = NULL
          ) }) %>%
        gather(stability_term, estimate, -trait, -marker, -b_std_error, -df) %>%
        ungroup() )


## Add marker information (position, etc.)
marker_stability_coef_info <- marker_stability_coef %>% 
  map(~left_join(rename(snp_info, marker = `rs#`), ., by = "marker"))

# Add the regression coefficients to the marker effects
S2_MET_marker_eff_pheno_fw <- list(S2_MET_marker_env_eff, marker_stability_coef_info) %>%
  pmap(~left_join(.x, .y, by = c("marker", "trait"))) %>%
  # Rearrange some columns and edit some variables
  map(~select(., trait, environment, env_effect = h, marker, chrom, alleles, pos, cM_pos, 
              mar_effect = effect, b_std_error, df, stability_term, estimate))


## Save the results
save_file <- file.path(result_dir, "S2MET_pheno_mar_eff_fw_results.RData")
save("S2_MET_marker_eff_pheno_fw", file = save_file)
