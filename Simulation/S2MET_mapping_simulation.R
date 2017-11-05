## S2MET Mapping Simulation
## 
## Run a simulation to test different association mapping models to detect QxE
## 

# Packages and directories
# List of packages
packages <- c("dplyr", "purrr", "tibble", "tidyr", "readr", "stringr", "readxl", "modelr", 
              "parallel", "purrrlyr", "rrBLUP", "pbr", "pbsim")

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

# Use only the tp_geno data
s2_imputed_mat_use <- s2_imputed_mat[c(tp_geno),]

# Round the genotype calls to be -1, 0, 1
panel_genos <- round(s2_imputed_mat_use) %>%
  ifelse(. < -1, -1) %>%
  ifelse(. > 1, 1)



### Create the genome and genetic model
gen_map <- snp_info %>% 
  select(1,3,5) %>% 
  data.frame(row.names = .$`rs#`) %>% 
  select(-1) %>% 
  table_to_map() %>% 
  jittermap()

genome <- sim_genome(map = gen_map)



### Set some simulation parameters
n_qtl <- c(30, 100)
h2 <- c(0.2, 0.5, 0.8)
n_iter <- 500

### Contstant parameters
n_env <- 5


### Simulate a genetic model
genome1 <- sim_gen_model(genome = genome, qtl.model = matrix(NA, nrow = n_qtl, ncol = 4),
                         add.dist = "geometric")

# Create a population
assoc_panel <- create_pop(genome = genome1, geno = panel_genos + 1)

# Phenotype
assoc_panel_pheno <- sim_phenoval(pop = assoc_panel, h2 = h2, n.env = n_env, n.rep = 1)

# Extract phenotypes
pheno <- assoc_panel_pheno$pheno_val$pheno_obs %>% 
  spread(trait, phenoval) %>% 
  select(-rep) %>%
  as.data.frame()

# Genotype
assoc_panel_geno <- genotype(genome = genome1, pop = assoc_panel) %>%
  t() %>%
  data.frame(marker = row.names(.), ., row.names = NULL, stringsAsFactors = FALSE,
             check.names = FALSE)

## Combine genotypes with snp information
geno <- right_join(x = snp_info, assoc_panel_geno, c("rs#" = "marker")) %>%
  rename(marker = `rs#`) %>%
  select(-alleles, -cM_pos) %>%
  as.data.frame()

### Association analysis
# K model
K_fit <- GWAS(pheno = pheno, geno = geno, fixed = "env", n.PC = 0, min.MAF = 0,
              plot = FALSE)

## Get the output




















