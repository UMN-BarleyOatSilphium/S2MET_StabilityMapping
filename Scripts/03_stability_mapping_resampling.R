## S2MET GWAS of Stability Coefficients
##
## Perform association analysis of the FW stability coefficients for regular phenotype
## FW, and for FW on environmental covariables
##
##

# Load FarmCPU
library("bigmemory")
require("biganalytics")
library("compiler") #this library is already installed in R
source("gapit_functions.R")
source("FarmCPU_functions.R")

library(snpStats)
library(GWASTools)


# Project and other directories
source(here::here("startup.R"))

## Association analysis
# Significance threshold
fdr_level <- 0.10


# Load data ---------------------------------------------------------------

# Load project data
load(file.path(data_dir, "project_data.RData"))

## Load stability data
load(file.path(result_dir, "stability_estimates.RData"))

# Reformat the genotype data for rrBLUP
geno_df <- as.data.frame(inner_join(select(snp_info, marker, chrom, pos),
                                    rownames_to_column(as.data.frame(t(geno_mat_mapping)), "marker")))

# Subset just the TP
geno_df_tp <- select(geno_df, marker, chrom, pos, all_of(tp_geno))


# Create a wide data.frame
pheno_df <- genotype_stability %>%
  select(trait, regression) %>%
  unnest(regression) %>%
  mutate(sd_d = sqrt(var_d)) %>%
  select(line_name, trait, g, b, sd_d) %>%
  gather(parameter, value, g, b, sd_d)

pheno_df_tp <- subset(pheno_df, line_name %in% tp_geno)

pheno_df_wide <- pheno_df %>%
  unite(trait.param, trait, parameter, sep = ".") %>%
  spread(trait.param, value) %>%
  as.data.frame()

pheno_df_wide_tp <- subset(pheno_df_wide, line_name %in% tp_geno)





# Convert the genotype matrix for FarmCPU
# Convert to [0,2] scale
geno_mat_mapping1 <- geno_mat_mapping + 1
geno_mat_mapping1 <- ifelse(geno_mat_mapping1 > 2, 2, ifelse(geno_mat_mapping1 < 0, 0, geno_mat_mapping1))
GD <- as.big.matrix(t(geno_mat_mapping1))
GM <- snp_info %>%
  select(marker, chrom, pos) %>%
  mutate(marker = factor(marker, levels = row.names(GD))) %>%
  arrange(marker) %>%
  as.data.frame()

# PCs for the K matrix
Ktp <- A.mat(X = geno_mat_mapping, min.MAF = 0, max.missing = 1)

Kpca <- broom::tidy(prcomp(x = Ktp)) %>%
  mutate(PC = paste0("PC", str_pad(string = PC, width = max(nchar(PC)), pad = "0"))) %>%
  spread(PC, value) %>%
  as.data.frame() %>%
  column_to_rownames("row") %>%
  as.matrix()

# Choose the best model

# Read in the filled version
univariate_gwas_farm_bestModel_qq <- read_csv(file = file.path(result_dir, "univariate_gwas_farm_bestModel_filled.csv"))



# Run GWAS for each stability sample --------------------------------------

# Create new working directory
farmcpu_dir <- file.path(result_dir, "FarmCPU_GWAS_Samples")
if (!dir.exists(farmcpu_dir)) {
  dir.create(farmcpu_dir)
}
setwd(dir = farmcpu_dir)


pheno_df_samples_tp <- genotype_stability_samples %>%
  select(trait, pEnv, rep, regression) %>%
  unnest(regression) %>%
  mutate(sd_d = sqrt(var_d)) %>%
  select(trait, line_name, pEnv, rep, g, b, sd_d) %>%
  gather(parameter, value, g, b, sd_d) %>%
  filter(line_name %in% tp_geno)

# Iterate over trait
univariate_gwas_sample_farm_out <- pheno_df_samples_tp %>%
  group_by(trait, parameter, pEnv, rep) %>%
  nest() %>%
  ungroup() %>%
  left_join(., select(univariate_gwas_farm_bestModel_qq, trait, parameter, nPC = bestModel_qq_nPC)) %>%
  mutate(gwas_out = list(NULL)) %>%
  arrange(trait, pEnv, parameter)

pb <- progress::progress_bar$new(total = sum(sapply(univariate_gwas_sample_farm_out$gwas_out, is.null)))

# Iterate over rows
for (i in which(sapply(univariate_gwas_sample_farm_out$gwas_out, is.null))) {
  
  # Create a new subdir
  trt <- univariate_gwas_sample_farm_out$trait[i]
  nEnvSam <- univariate_gwas_sample_farm_out$pEnv[i]
  param <- univariate_gwas_sample_farm_out$parameter[i]
  rep <- univariate_gwas_sample_farm_out$rep[i]
  
  # Rename the trait
  subdir <- paste0(trt, "_pEnv", sub(pattern = "\\.", replacement = "", x = nEnvSam), "_", param)
  new_trait_name <- paste0(subdir, "_", rep)
  
  # Create the subdir
  full_subdir <- file.path(farmcpu_dir, subdir)
  if (!dir.exists(full_subdir)) {
    dir.create(full_subdir)
  }
  
  # Change working directory
  setwd(full_subdir)
  
  dat <- univariate_gwas_sample_farm_out$data[[i]] %>%
    droplevels() %>%
    arrange(line_name)
  
  Y <- as.data.frame(dat)
  
  # Number of PC covariates
  n_PC <- univariate_gwas_sample_farm_out$nPC[i]
  if (n_PC == 0) {
    PC_covar <- NULL
  } else {
    PC_covar <- Kpca[levels(Y$line_name), seq_len(n_PC), drop = FALSE]
  }
  
  names(Y)[2] <- new_trait_name
  
  # Try running
  try_fcpu <- try({
    out <- capture.output({
      farmCPU_out <- FarmCPU(Y = Y, GD = GD[,levels(Y$line_name)], GM = GM, CV = PC_covar, maxLoop = 10, MAF.calculate = TRUE,
                             maf.threshold = 0)
    })
  }, silent = TRUE)
  
  # If it error-ed, skip
  if (inherits(try_fcpu, "try-error")) {
    next
    
  } else {
    # Calculate the FDR threshold
    fdr_p_thresh <- sommer:::fdr(p = farmCPU_out$GWAS$P.value, fdr.level = fdr_level)
    farmCPU_out$fdr_p_thresh <- fdr_p_thresh$fdr.10
    
    univariate_gwas_sample_farm_out$gwas_out[[i]] <- farmCPU_out
    
  }
  
  pb$tick()
  
}


univariate_gwas_sample_farm <- univariate_gwas_sample_farm_out %>%
  select(-data)

# Save everything
save("univariate_gwas_sample_farm", file = file.path(result_dir, "stability_samples_gwas_out.RData"))






