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
pheno_dir <- file.path(alt_proj_dir, "Phenotype_Data/")
env_var_dir <- file.path(alt_proj_dir, "Environmental_Variables")

geno_dir <-  "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/Data/GBS_Genos"
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

# Extract the tp and vp from the G matrix
s2_imputed_mat_use <- s2_imputed_mat[c(tp_geno, vp_geno),]
# s2_imputed_mat_use <- s2_imputed_mat[tp_geno,]

# Remove the environments in which the vp was only observed
S2_MET_BLUEs_use <- S2_MET_BLUEs %>%
  filter(!grepl(pattern = "BZI|HTM", x = environment)) %>%
  group_by(trait, environment) %>%
  filter(sum(line_name %in% tp) > 1,
         line_name %in% c(tp_geno, vp_geno)) %>%
  ungroup() %>%
  mutate(line_name = as.factor(line_name),
         environment = as.factor(environment))

## Calculate marker effects x environment as in Lopez-Cruz 2015
# 
# # Test on three environments
# test_blues <- S2_MET_BLUEs_use %>%
#   filter(trait == "GrainYield", environment %in% c("BCW16", "STP16", "CRM16"))
# 
# # Training data
# train_blues <- test_blues %>%
#   filter(line_name %in% tp_geno) %>%
#   droplevels()

# Marker matrix
M <- s2_imputed_mat_use

mar_names <- colnames(M)


## Group by trait and estimate the base marker effect and enivornmnet-specific deviation
mar_effect_by_env <- S2_MET_BLUEs_use %>%
  group_by(trait) %>%
  do({
    
    # Convert the . to df
    df <- .
    # Create a model frame with lines and environments
    mf <- model.frame(value ~ line_name + environment, data = df, drop.unused.levels = TRUE)
    env_names <- levels(mf$environment)

    # Response vector
    y <- model.response(mf)

    # Model matrix for the fixed effects
    X <- model.matrix(~ -1 + environment, data = mf)

    # Model matrix of the line_name incidence
    Z_line <- model.matrix(~ -1 + line_name, mf)
    # Model matrix of the marker random effects
    
    ## First the main effect of each marker
    Z0 <- Z_line %*% M
    K0 <- Diagonal(ncol(Z0))
    
    # Random deviation of each marker in each environment
    Z1 <- mf %>% 
      split(.$environment) %>% 
      map(~ model.matrix(~ -1 + line_name, .)) %>% 
      map(~ . %*% M)
    Z1 <- .bdiag(Z1)
    
    K1 <- Diagonal(ncol(Z1))
    
    fit <- EMMREML::emmremlMultiKernel(y = y, X = X, Zlist = list(Z0, Z1), Klist = list(K0, K1))

    # Rename the random effect vector
    fit$uhat %>%
      setNames(c(mar_names, paste(rep(mar_names, length(env_names)), 
                                  rep(env_names, each = length(mar_names)), sep = ":"))) %>%
      data.frame(marker = names(.), effect = ., row.names = NULL, stringsAsFactors = FALSE) %>%
      separate(marker, c("marker", "environment"), sep = ":", fill = "right")
    
  })

## Save this
save_file <- file.path(result_dir, "S2MET_marker_eff_by_env.RData")
save("mar_effect_by_env", file = save_file)


## Load the marker effect by environment results
load(file.path(result_dir, "S2MET_marker_eff_by_env.RData"))
# Loa the FW regression results
load(file.path(result_dir, "S2MET_fw_regression_results.RData"))


## Combine the marker effects by environment with the environmental effect from
## FW regression
S2_MET_marker_env <- left_join(filter(mar_effect_by_env, !is.na(environment)), 
                               distinct(S2_MET_BLUEs_fw, trait, environment, h))

# Run FW regression of marker effects per environment
S2_MET_marker_fw <- S2_MET_marker_env %>% 
  group_by(trait, marker) %>% 
  # Need to adjust the scale so they are the same - remember for genotypes and environments
  # they are already on the same scale 
  # mutate_at(vars(effect, h), scale) %>% 
  do({
    df <- .
    # Fit the model
    fw_fit <- lm(effect ~ h, data = df)
    df_resid <- df.residual(fw_fit)
    # Grab the coefficients and standard errors
    fw_coef <- tidy(fw_fit) %>% 
      subset(term == "h", term:std.error) %>%
      # Calculate the resid MS of the regression (E and R 1966 stability)
      add_row(term = "delta", estimate = sum(resid(fw_fit)^2) / df_resid) %>%
      mutate(df = df_resid)
    
    # Return a data.frame with the model fit, the coefficients, and the residuals
    data_frame(coef = list(fw_coef), resid = list(resid(fw_fit)))
  })

# Calculate the pooled error
S2_MET_marker_fw1 <- S2_MET_marker_fw %>% 
  group_by(trait) %>% 
  # Sum of squared residuals / residual df
  mutate(pooled_df = (length(unlist(resid)) - (2 * n())),
         pooled_error = sum(unlist(resid)^2) / pooled_df)

# Run significance tests for the "b" term and the "delta" term
S2_MET_marker_fw_sig <- S2_MET_marker_fw1 %>% 
  unnest(coef) %>%
  mutate(term = if_else(term == "h", "b", term),
         std.error = if_else(term == "b", std.error, pooled_error),
         statistic = estimate / std.error,
         # The test for stability is a t test for the regression coefficient and
         # an F test for the deviations from the predictions
         stable = if_else(term == "b", 1 - pt(q = statistic, df = df, lower.tail = FALSE),
                            1 - pf(q = statistic, df1 = df, df2 = pooled_df, lower.tail = FALSE)),
         sensitive = 1 - stable) %>%
  select(-pooled_df, -pooled_error, -df) %>%
  gather(test, p_value, -trait:-statistic) %>%
  ungroup()
    

## Save the results
save_file <- file.path(result_dir, "S2MET_marker_eff_fw_results.RData")
save("S2_MET_marker_fw_sig", file = save_file)
