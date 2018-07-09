## Source script for the S2MET_Mapping project
## 
## 

## Load commonly used libraries
library(tidyverse)
library(readxl)
library(rrBLUP)

# Load personal libraries
library(neyhart)


## The Google Drive directory
# gdrive_dir <- "C:/Users/Jeff/GoogleDrive"
gdrive_dir <- "/Users/neyha001/Google Drive/"

# Project and other directories
proj_dir <- repo_dir
alt_proj_dir <- file.path(gdrive_dir, "BarleyLab/Projects/S2MET")

# Geno, pheno, and enviro data
geno_dir <-  file.path(gdrive_dir, "BarleyLab/Projects/Genomics/Genotypic_Data/GBS_Genotype_Data/")
# geno_dir <- 
bopa_geno_dir <- file.path(gdrive_dir, "/BarleyLab/Projects/Genomics/Genotypic_Data/BOPA_Genotype_Data/")
# bopa_geno_dir <- 
pheno_dir <- file.path(alt_proj_dir, "Data/")

# Character vector for replacing stability/mean variable symbols with full name
coef_replace <- c("b" = "Linear Stability",
                  "log_delta" = "Non-Linear Stability",
                  "g" = "Genotype Mean")


# Color scheme for manhattan plot
color <- c(set_names(umn_palette(n = 4)[3:4], "Bl", "Or"), "B" = "black", "G" = "grey75")

# Placeholder
alpha <- 0

### Plotting modifiers
# Manhattan plot
g_mod_man <- list(
  geom_point(),
  # geom_hline(aes(yintercept = neg_log10_fdr10, lty = "FDR 10%")),
  scale_color_manual(values = color, guide = FALSE),
  ylab(expression(-log[10](italic(q)))),
  xlab("Position (Mbp)"),
  theme_bw(),
  theme_manhattan(),
  theme(panel.border = element_blank())
)


######
# MSI Source starts here
######

# Source the project functions
source(file.path(proj_dir, "source_functions.R"))

# Other directories
script_dir <- file.path(proj_dir, "Scripts/")
fig_dir <- file.path(proj_dir, "Figures/")
data_dir <- file.path(proj_dir, "Data")

map_dir <- file.path(script_dir, "GWAS")
result_dir <- file.path(proj_dir, "Results")

# Load the phenotypic data
load(file.path(pheno_dir, "S2_MET_BLUEs.RData"))
# Load the genotypic data
load(file.path(bopa_geno_dir, "S2TP_multi_genos.RData"))

# Load the genotype data for the TP and VP
load(file.path(geno_dir, "S2_genos_hmp.RData"))
load(file.path(geno_dir, "S2_genos_mat.RData"))



# Load an entry file
entry_list <- read_excel(file.path(data_dir, "S2MET_project_entries.xlsx"))


# Grab the entry names that are not checks
tp <- subset(entry_list, Class == "S2TP", Line, drop = TRUE)

# Grab the vp
vp <- subset(entry_list, Class == "S2C1R", Line, drop = TRUE)
  

# Find the tp lines in the genotype matrix
tp_geno <- intersect(tp, row.names(S2TP_imputed_multi_genos_mat))

# Find the tp and vp lines in the other genotype matrix
tp_geno1 <- intersect(tp, row.names(s2_imputed_mat))
vp_geno1 <- intersect(vp, row.names(s2_imputed_mat))




# Define the checks
checks <- subset(entry_list, Class == "Check", Line, drop = TRUE)

entries <- entry_list$Line

# Filter the BLUEs for the TP
S2_MET_BLUEs_tp <- S2_MET_BLUEs %>% 
  filter(line_name %in% tp,
         trait %in% c("GrainYield", "HeadingDate", "PlantHeight"))

# Filter the BLUEs for the TP+VP
S2_MET_BLUEs_tpvp <- S2_MET_BLUEs %>% 
  filter(line_name %in% c(tp, vp),
         trait %in% c("GrainYield", "HeadingDate", "PlantHeight"))


# Read in the trial metadata
trial_info <- read_csv(file = file.path(file.path(alt_proj_dir, "Data/"), "trial_metadata.csv"))


## Create a snp information df
snp_info <- S2TP_imputed_multi_genos_hmp %>% 
  select(marker = 1, chrom:cM_pos)

# Same for GBS snps
snp_info1 <- s2_discrete_genos %>%
  select(marker = 1, chrom:cM_pos)



