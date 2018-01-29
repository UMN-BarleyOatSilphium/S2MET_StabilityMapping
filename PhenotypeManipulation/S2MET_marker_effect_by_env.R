## S2MET Mapping
## Calculate marker effects in each environment
## 
## This script will calculate marker effects across environment using a G model, similar to GWAS.
## 

# Capture the trait to use (argument 1)
args <- commandArgs(trailingOnly = T)
tr <- args[1]

# List of packages to load
# List of packages
packages <- c("dplyr", "purrr", "tibble", "tidyr", "readr", "stringr", "readxl", "modelr", 
              "parallel", "purrrlyr", "rrBLUP", "ggplot2", "broom", "Matrix", "pbr", "lme4qtl")

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