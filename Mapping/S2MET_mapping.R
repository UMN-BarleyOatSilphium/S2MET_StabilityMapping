## Perform GWAS using the S2MET data
## 
## GWAS will be performed using the 1) phenotype BLUEs with environments as a fixed
## effect (i.e. phenotypic mean mapping), 2) all phenotypes to detect QTLxE interaction,
## 3) the main effect and stability coefficient from FW regression


# List of packages to load
# List of packages
packages <- c("tidyverse", "stringr", "readxl", "modelr", "parallel", "purrrlyr",
              "rrBLUP", "sommer")

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
geno_dir <-  "C:/Users/Jeff/Google Drive/Barley Lab/Projects/Genomic Selection/Genotypic Data/GBS Genotype Data/"
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


### Setup

# Calculate the A matrix
A <- A.mat(X = s2_imputed_mat_use, min.MAF = 0, max.missing = 1)

## Visualize pop structure
A_prcomp <- prcomp(A)

# Extract the proportion explained by each PC
lambda <- round(A_prcomp$sdev / sum(A_prcomp$sdev), 3) * 100

# Convert to data.frame
A_prcomp_df <- A_prcomp$rotation %>% 
  as.data.frame() %>% 
  rownames_to_column("line_name")

A_prcomp_df1 <- left_join(A_prcomp_df, entry_list, c("line_name" = "Line")) %>% 
  dplyr::select(line_name, program = Program, parent = `Parent?`, starts_with("PC")) %>%
  mutate(program = ifelse(str_detect(line_name, "^2MS"), "M2", program))

# Plot
(g_pop_str <- A_prcomp_df1 %>% 
  ggplot(aes(x = PC1, y = PC2, col = program)) + 
  geom_point(size = 2) + 
  ylab(str_c("PC2 (", lambda[2], "%)")) +
  xlab(str_c("PC1 (", lambda[1], "%)")) +
  scale_color_discrete(guide = guide_legend(title = "Breeding\nProgram")) +
  labs(
    title = "Population Structure of the S2MET Genotypes",
    subtitle = "The 'M2' program includes genotypes derived from crosses made among genotypes from the other programs."
  ))
  




### GWAS with Phenotypic Means

# Prepare phenotypic data
# Filter out the environments in which only the C1 population was grown
phenos_use <- S2_MET_BLUEs %>%
  group_by(trait, environment) %>% 
  filter(n() > 50) %>% 
  ungroup() %>%
  subset(select = c(line_name, environment, trait, value)) %>%
  spread(trait, value) %>%
  filter(line_name %in% c(tp_geno, vp_geno)) %>%
  as.data.frame()

# Prepare genotypic data
genos_use <- s2_imputed_genos %>%
  dplyr::select(marker = `rs#`, chrom, pos, which(names(.) %in% c(tp_geno, vp_geno))) %>%
  as.data.frame()

# Randomly sample 5 environments
sample_env <- sample(unique(phenos_use$environment), 5)

## Test with PlantHeight
phenos_test <- subset(phenos_use, environment %in% sample_env,
                      select = c(line_name, environment, PlantHeight))


# Run GWAS
gwas_main_out <- rrBLUP::GWAS(pheno = phenos_test, geno = genos_use, fixed = "environment",
                      n.PC = 2, min.MAF = 0, P3D = FALSE, K = A)

## Create model matrices
mf <- model.frame(formula = PlantHeight ~ line_name + environment, data = phenos_test, )
y <- model.response(mf)

# Fixed effects of environments
X <- model.matrix(~ environment, data = mf)
# Random effect of genotype
Z <- gws::ranef_model_matrix(~ g(line_name), data = mf, vcov = list(line_name = A))

# Fixed effect of markers
W <- s2_imputed_mat_use


# Test sommer
gwas_main_out_somm <- GWAS(Y = y, X = X, Z = Z, W = W, min.MAF = 0, gwas.plots = FALSE)

# Save
save_file <- file.path(map_dir, "S2MET_gwas_main_results.RData")
save("gwas_main_out", file = save_file)



