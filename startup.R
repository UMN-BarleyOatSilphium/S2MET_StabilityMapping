<<<<<<< HEAD
## Source script for the S2MET_Mapping project
## 
## 

## Load commonly used libraries
library(tidyverse)
library(readxl)
library(rrBLUP)

library(neyhart)

# Project and other directories
proj_dir <- here::here()

# Other directories
fig_dir <- file.path(proj_dir, "Figures")
data_dir <- file.path(proj_dir, "Data")
result_dir <- file.path(proj_dir, "Results")





# Character vector for replacing stability/mean variable symbols with full name
# coef_replace <- c("b" = "Linear Stability",
#                   "log_delta" = "Non-Linear Stability",
#                   "g" = "Genotype Mean")

coef_replace <- c("g" = "Genotype Mean", "b" = "Linear Stability", "sd_d" = "Non-Linear Stability")
coef_replace1 <- c("g" = "Mean", "b" = "Slope", "sd_d" = "Deviation")

# Color scheme for manhattan plot
color <- c(set_names(neyhart_palette(name = "umn1", n = 4)[3:4], "Bl", "Or"), "B" = "black", "G" = "grey75")

# Significance level
alpha <- 0.05

# Traits of interest
traits <- c("GrainYield", "HeadingDate", "PlantHeight", "TestWeight", "GrainProtein")

# Trait units
trait_units <- c("kg ha^-1", "days", "cm", "g L^-1", "'%'")
trait_units <- setNames(trait_units, traits)

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

# Source the project functions'
source(file.path(proj_dir, "functions.R"))


=======
## Source script for the S2MET_Mapping project
## 
## 

## Load commonly used libraries
library(tidyverse)
library(readxl)
library(rrBLUP)

library(neyhart)

# Project and other directories
proj_dir <- here::here()

root <- neyhart::find_dir("SideProjects")

# Geno, pheno, and enviro data
geno_dir <-  file.path(root, "ProjectData/GenotypicData/")
pheno_dir <- file.path(root, "ProjectData/PhenotypicData/")
meta_dir <- pheno_dir
enviro_dir <- file.path(root, "ProjectData/EnvironmentalData/")

# Other directories
fig_dir <- file.path(proj_dir, "Figures")
data_dir <- file.path(proj_dir, "Data")
result_dir <- file.path(proj_dir, "Results")





# Character vector for replacing stability/mean variable symbols with full name
# coef_replace <- c("b" = "Linear Stability",
#                   "log_delta" = "Non-Linear Stability",
#                   "g" = "Genotype Mean")

coef_replace <- c("g" = "Genotype Mean", "b" = "Linear Stability", "sd_d" = "Non-Linear Stability")


# Color scheme for manhattan plot
color <- c(set_names(neyhart_palette(name = "umn1", n = 4)[3:4], "Bl", "Or"), "B" = "black", "G" = "grey75")

# Significance level
alpha <- 0.05

# Traits of interest
traits <- c("GrainYield", "HeadingDate", "PlantHeight", "TestWeight", "GrainProtein")

# Trait units
trait_units <- c("kg ha^-1", "days", "cm", "g L^-1", "'%'")
trait_units <- setNames(trait_units, traits)

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

# Source the project functions'
source(file.path(proj_dir, "functions.R"))


>>>>>>> c042eb08f49fab9587e31a12e4f60ec19baa31ee
