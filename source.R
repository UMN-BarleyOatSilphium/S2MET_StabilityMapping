## Source script for the S2MET_Mapping project
## 
## 

# Project and other directories
proj_dir <- repo_dir
alt_proj_dir <- "C:/Users/Jeff/Google Drive/Barley Lab/Projects/S2MET"

# Geno, pheno, and enviro data
geno_dir <-  "C:/Users/Jeff/Google Drive/Barley Lab/Projects/Genomics/Genotypic_Data/GBS_Genotype_Data/"
# geno_dir <- 
bopa_geno_dir <- "C:/Users/Jeff/Google Drive/Barley Lab/Projects/Genomics/Genotypic_Data/BOPA_Genotype_Data/"
# bopa_geno_dir <- 
pheno_dir <- file.path(alt_proj_dir, "Phenotype_Data/")

# Other directories
script_dir <- file.path(proj_dir, "Scripts/")
analysis_dir <- file.path(script_dir, "Analysis")
fig_dir <- file.path(script_dir, "Figures/")

map_dir <- file.path(analysis_dir, "GWAS")
result_dir <- file.path(proj_dir, "Results")

# Load the phenotypic data
load(file.path(pheno_dir, "S2_MET_BLUEs.RData"))
# Load the genotypic data
load(file.path(bopa_geno_dir, "S2TP_multi_genos.RData"))

# Load an entry file
entry_list <- read_excel(file.path(proj_dir, "S2MET_project_entries.xlsx"))


# Grab the entry names that are not checks
tp <- entry_list %>% 
  filter(Class == "S2TP") %>% 
  pull(Line)


# Find the tp and vp that are genotypes
tp_geno <- intersect(tp, row.names(S2TP_imputed_multi_genos_mat))

# Define the checks
checks <- entry_list %>% 
  filter(Class == "Check") %>% 
  pull(Line)

entries <- entry_list %>% 
  pull(Line)

# Filter environments for those in which the TP was observed
S2_MET_BLUEs_use <- S2_MET_BLUEs %>% 
  filter(line_name %in% tp_geno,
         trait %in% c("GrainYield", "HeadingDate", "PlantHeight"))

# Read in the trial metadata
trial_info <- read_csv(file = file.path(pheno_dir, "trial_metadata.csv"))

# Character vector for replacing stability/mean variable symbols with full name
coef_replace <- c("b" = "Linear Stability", "log_delta" = "Non-Linear Stability",
                  "g" = "Genotype Mean")
