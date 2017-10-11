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
# load(file.path(pheno_dir, "S2_MET_BLUEs.RData"))
# Load the genotypic data
load(file.path(geno_dir, "S2_genos_mat.RData"))
# load(file.path(geno_dir, "S2_genos_hmp.RData"))

# Load environmental data
# load(file.path(env_var_dir, "environmental_data_compiled.RData"))

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


# Load the FW results
load(file.path(result_dir, "S2MET_fw_regression_results.RData"))

# Format the phenotypic data and only select the TP
S2_MET_BLUEs_fw_crossv <- S2_MET_BLUEs_fw %>% 
  filter(line_name %in% c(tp_geno)) %>%
  select(line_name, trait, g, b) %>%
  distinct() %>%
  gather(coef, value, -line_name, -trait)


# Set the number of CV iterations and the proportion to use as training
n_cv_iter <- 1000
p_train <- 0.70


# Create training-test combinations of the modeling data
cv_parts <- S2_MET_BLUEs_fw_crossv %>% 
  group_by(trait, coef) %>% 
  do(cv_parts = crossv_mc(data = ., n = n_cv_iter, test = 1 - p_train)) %>%
  ungroup()

# Create a function to output the prediction accuracy
cv_pred_acc <- function(geno_train, geno_test, pheno_train, pheno_test) {
  
  # Model frame
  mf <- model.frame(value ~ line_name, data = pheno_train)
  
  # y vector
  y <- model.response(mf)
  # Fit the model
  fit <- mixed.solve(y = y, Z = geno_train, method = "REML")
  
  # Extract the marker effects and calculate the PGVs
  PGV <- geno_test %*% fit$u
  
  # Calculate accuracy
  mf_test <- model.frame(value ~ line_name, pheno_test)
  
  y_test <- model.response(mf_test)
  as.numeric(cor(PGV, y_test))
  
}

# Group by trait and coefficient and then run the cross-validation
crossv_gbs_out <- cv_parts %>% 
  group_by(trait, coef) %>% 
  do({
    
    mcs <- .$cv_parts[[1]]
    
    cv_accs <- apply(X = mcs, MARGIN = 1, FUN = function(i) {
      
      # Designate phenos for training and testing
      pheno_train <- as.data.frame(i$train)
      pheno_test <- as.data.frame(i$test)
      
      # Designate genos for training and testing
      geno_train <- s2_imputed_mat_use[pheno_train$line_name,]
      geno_test <- s2_imputed_mat_use[pheno_test$line_name,]
      
      # Predict
      cv_pred_acc(geno_train = geno_train, geno_test = geno_test,
                  pheno_train = pheno_train, pheno_test = pheno_test)
      
    }) 
    
    # Return a dataframe
    data_frame(acc = mean(cv_accs), acc_sd = sd(cv_accs))
    
    })

crossv_gbs_out <- crossv_gbs_out %>% 
  mutate(crossv_iter = n_cv_iter, p_train = p_train)

# Save this
save_file <- file.path(result_dir, "S2MET_stability_crossv_results.RData")
save("crossv_gbs_out", file = save_file)















