## S2MET GWAS of Stability Coefficients - resampling
## 
## Perform association analysis of the FW stability coefficients for regular phenotype
## FW, and for FW on environmental covariables
## 
## Author: Jeff Neyhart
## Last updated: May 25, 2018
## 


# Run the source script - MSI
repo_dir <- "/panfs/roc/groups/6/smithkp/neyha001/QTLMapping/S2MET_Mapping/"
source(file.path(repo_dir, "source_MSI.R"))

# Run the source script - local
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))




## Load packages
packages <- c("lme4")
invisible(lapply(packages, library, character.only = TRUE, lib.loc = package_dir))