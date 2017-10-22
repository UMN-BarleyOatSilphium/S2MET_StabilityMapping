## S2MET Mapping
## Identify and map marker effect stability
## 
## This script will calculate marker effects across environment, test for significant
## marker effect stability, then attempt to map markers will stable or sensitive effects
## 

# List of packages to load
# List of packages
packages <- c("dplyr", "purrr", "tibble", "tidyr", "readr", "stringr", "readxl", "modelr", 
              "parallel", "purrrlyr", "rrBLUP", "FW", "ggplot2", "broom")

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
load(file.path(geno_dir, "S2_genos_hmp.RData"))

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


## Load the FW regression results
load(file.path(result_dir, "S2MET_fw_regression_results.RData"))


# Remove the environments in which the vp was only observed
S2_MET_BLUEs_use <- S2_MET_BLUEs %>%
  group_by(trait, environment) %>%
  filter(sum(line_name %in% tp_geno) > 1,
         line_name %in% c(tp_geno, vp_geno)) %>%
  ungroup() %>%
  mutate(line_name = as.factor(line_name))

S2_MET_BLUEs_fw <- S2_MET_BLUEs_use %>% 
  group_by(trait) %>%
  do({
    FW(y = .$value, VAR = .$line_name, ENV = .$environment, method = "OLS")$h %>%
      as.data.frame() %>% rownames_to_column("environment") %>% rename(h = V1)
  })


# Compound symmetry + diagonal model for marker effects
# First return the marker matrix for each environment
S2_MET_BLUEs_use %>%
  group_by(trait, environment) %>%
  do({
    
    df <- .

    # Create a model.frame
    mf <- model.frame(value ~ line_name, data = df)
    # Get the model response
    y <- model.response(mf)
    # Get the matrix of line_name observations
    z <- model.matrix(~ -1 + line_name, mf)

## Predict marker effects
# For each environment and trait, predict the marker effects
marker_effects_by_env <- S2_MET_BLUEs_use %>% 
  # filter(line_name %in% tp_geno) %>% droplevels() %>%
  group_by(trait, environment) %>%
  do({
    
    df <- .
    
    # Create a model.frame
    mf <- model.frame(value ~ line_name, data = df)
    # Get the model response
    y <- model.response(mf)
    # Get the matrix of line_name observations
    z <- model.matrix(~ -1 + line_name, mf)
    
    
    
    # Extract observations from the marker matrix
    Z <- z %*% s2_imputed_mat_use
    # Z <- z %*% s2tp_bopa_mat_use
    
    
    # Fit the model
    fit <- mixed.solve(y = y, Z = Z, method = "REML", SE = FALSE, return.Hinv = FALSE)
    
    # Extract and return the marker effects
    fit$u %>% 
      data_frame(marker = names(.), effect = .)
    
  })
    

# Add the gradients to the marker effect environment data frame
marker_effects_by_env1 <- marker_effects_by_env %>% 
  left_join(., distinct(S2_MET_BLUEs_fw, trait, environment, h), by = c("trait", "environment")) %>%
  # Add environmental variables
  left_join(., spread(one_year_env_df, variable, value)) %>%
  gather(gradient, value, -trait:-effect)


# Now fit the FW regression model for each marker using h or environmental gradients as the predictor
marker_fw_fit <- marker_effects_by_env1 %>%
  group_by(trait, marker, gradient) %>% 
  do({
    fw_fit <- lm(effect ~ value, data = .)
    data_frame(tidy = list(tidy(fw_fit)), df = df.residual(fw_fit))
  })
    


# Extract the coefficient for h and perform a signficance test for greater
# than zero or less than zero
marker_fw_coef <- marker_fw_fit %>% 
  unnest() %>% 
  filter(term == "h") %>%
  mutate(stable = pt(statistic, df = df, lower.tail = TRUE),
         sensitive = pt(statistic, df = df, lower.tail = FALSE)) %>%
  gather(test, p_value, -trait:-p.value)


## Map
# Adjust the p-values for multiple testing and convert to -log10
marker_fw_adj <- marker_fw_coef %>% 
  group_by(trait, test) %>% 
  mutate(p_val_adj = p.adjust(p_value, method = "fdr", n = n()), 
         neg_log_p = -log10(p_val_adj),
         fdr05 = -log10(0.05),
         fdr10 = -log10(0.10))


# Add marker position data
marker_fw_adj1 <- s2_imputed_genos %>% 
  select(marker = `rs#`, alleles:pos) %>%
  mutate(chrom = as.factor(chrom)) %>%
  right_join(., marker_fw_adj, by = "marker") %>%
  arrange(trait, chrom, pos)


## plot
marker_fw_adj1 %>%
  ggplot(aes(x = pos, y = neg_log_p, col = chrom)) +
  geom_point() +
  facet_grid(trait + test ~ chrom, scales = "free") +
  geom_hline(aes(yintercept = fdr05, lty = "FDR 05%")) +
  geom_hline(aes(yintercept = fdr10, lty = "FDR 10%"))
  

