## S2MET GWAS of Stability Coefficients
##
## Perform association analysis of the FW stability coefficients for regular phenotype
## FW, and for FW on environmental covariables
##
##

library(sommer) # For univariate and multi-variate GWAS

# Load FarmCPU
library("bigmemory")
library("biganalytics")
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
  filter(trait %in% traits) %>%
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



# Run univariate GWAS -----------------------------------------------------

# # Iterate over trait and parameter
# univariate_gwas_out <- pheno_df_tp %>%
#   group_by(trait, parameter) %>%
#   nest() %>%
#   ungroup() %>%
#   crossing(., nPC = 1:3) %>%
#   mutate(gwas_out = list(NULL))
# 
# for (i in seq_len(nrow(univariate_gwas_out))) {
#   
#   dat <- univariate_gwas_out$data[[i]] %>%
#     droplevels()
#   
#   n_PC <- univariate_gwas_out$nPC[i]
#   
#   # fit
#   gwas_out <- sommer::GWAS(fixed = value ~ 1, random = ~ vsr(line_name, Gu = K1), n.PC = n_PC,
#                            data = dat, M = geno_mat_mapping_tp, gTerm = "u:line_name")
#   
#   # Calculate the FDR threshold
#   fdr_p_thresh <- sommer:::fdr(p = gwas_out$scores[,1], fdr.level = fdr_level)$fdr.10
#   
#   # Merge scores with SNP.info
#   scores <- gwas_out$scores %>%
#     as.data.frame() %>%
#     rownames_to_column("marker") %>%
#     rename(score = value) %>%
#     mutate(score = ifelse(is.infinite(score), as.numeric(NA), score))
#   
#   # Return the results
#   univariate_gwas_out$gwas_out[[i]] <- tibble(scores = list(scores), fdr_p_thresh)
#   
# }
# 
# univariate_gwas <- univariate_gwas_out %>%
#   select(-data) %>%
#   unnest(gwas_out)
# 
# 
# # calculate variance inflation factor for each model
# univariate_gwas_vif <- univariate_gwas %>%
#   mutate(p_value_infl = map_dbl(scores, ~QCEWAS::P_lambda(p = 10^-.x$score))) %>%
#   # Assign the model with inflation factor closes to 1 as the best model
#   group_by(trait, parameter) %>%
#   mutate(best_model = nPC[which.min(abs(1 - p_value_infl))],
#          best_model = ifelse(nPC == best_model, "*", "")) %>%
#   ungroup()
# 
# 
# # nPC = 3 is most often the best model; use that.
# 
# n_PC_use <- 3



## Univariate GWAS using FarmCPU

# Set working directory
farmcpu_dir <- file.path(result_dir, "FarmCPU_GWAS")
setwd(dir = farmcpu_dir)


# Iterate over trait
univariate_gwas_farm_out <- pheno_df_tp %>%
  group_by(trait, parameter) %>%
  nest() %>%
  ungroup() %>%
  crossing(., nPC = 0:3) %>%
  mutate(gwas_out = list(NULL))

pb <- progress::progress_bar$new(total = sum(sapply(univariate_gwas_farm_out$gwas_out, is.null)))

# Iterate over rows
for (i in which(sapply(univariate_gwas_farm_out$gwas_out, is.null))) {
  
  dat <- univariate_gwas_farm_out$data[[i]] %>%
    droplevels() %>%
    arrange(line_name)
  
  Y <- as.data.frame(dat)
  
  # Number of PC covariates
  n_PC <- univariate_gwas_farm_out$nPC[i]
  if (n_PC == 0) {
    PC_covar <- NULL
  } else {
    PC_covar <- Kpca[levels(Y$line_name), seq_len(n_PC), drop = FALSE]
  }
  
  # Rename the trait
  new_trait_name <- paste0(univariate_gwas_farm_out$trait[i], "_", univariate_gwas_farm_out$parameter[i], "_nPC", n_PC)
  names(Y)[2] <- new_trait_name
  
  # Try
  try_farm <- try({
    output <- capture.output({
      farmCPU_out <- FarmCPU(Y = Y, 
                             GD = GD[,levels(Y$line_name)], 
                             GM = GM, 
                             CV = PC_covar, 
                             maxLoop = 10, 
                             method.bin = "optimum",
                             MAF.calculate = TRUE,
                             maf.threshold = 0)
    })
  }, silent = TRUE)
  
  # Control flow if try_farm failed
  if (class(try_farm) != "try-error") {
    # Calculate the FDR threshold
    fdr_p_thresh <- sommer:::fdr(p = farmCPU_out$GWAS$P.value, fdr.level = fdr_level)
    farmCPU_out$fdr_p_thresh <- fdr_p_thresh$fdr.10
    
    # Get the sequence of QTN used as kinship covariates
    idx <- max(which(grepl(pattern = "seqQTN", x = output))) + 1
    seqQTN <- eval(parse(text = paste0("c(", str_replace_all(string = str_trim(str_remove(string = str_trim(output[idx]), pattern = "\\[1\\]")), pattern = "[ ]{1,}", replacement = ", "), ")")))
    
    kin_markers <- row.names(GD)[seqQTN]
    
    farmCPU_out$kinship.markers <- kin_markers
    
    univariate_gwas_farm_out$gwas_out[[i]] <- list(farmCPU_out)
    
  }

  pb$tick()
  
}


univariate_gwas_farm <- univariate_gwas_farm_out %>%
  select(-data) %>%
  # Remove null results
  filter(!sapply(gwas_out, is.null)) %>%
  # Unnest
  mutate(gwas_out = map(gwas_out, 1)) %>%
  mutate(gwas_out = map(gwas_out, ~as_tibble(map(.x, list)) %>% 
                          select_if(~!sapply(., function(x) any(is.null(x)))) %>% unnest(fdr_p_thresh)) ) %>%
  unnest(gwas_out) %>%
  rename(scores = GWAS, pred = Pred)


## Plot QQ plots
univariate_gwas_farm_qqplot_list <- univariate_gwas_farm %>%
  split(list(.$trait, .$parameter))

# Iterate over this list
for (i in seq_along(univariate_gwas_farm_qqplot_list)) {
  df <- univariate_gwas_farm_qqplot_list[[i]]
  trt <- unique(df$trait)
  param <- unique(df$parameter)
  filename <- paste0("GWAS_QQ_", trt, "_", param, ".png")
  
  ncol <- 2
  nrow <- ceiling(nrow(df) / ncol)
  
  png(filename = file.path(fig_dir, filename), width = 8, height = 4 * nrow, units = "in", res = 300)
  par(mfrow = c(nrow, ncol))
  # Iterate over data to plot
  for (j in seq_len(nrow(df))) {
    col <- df$nPC[j] + 1
    qqPlot(pval = df$scores[[j]]$P.value, col = col, main = paste0(trt, "; param: ", param, "; nPC: ", df$nPC[j]))
  }
  dev.off()
  
}


# calculate variance inflation factor for each model
univariate_gwas_farm_vif <- univariate_gwas_farm %>%
  mutate(p_value_infl = map_dbl(scores, ~QCEWAS::P_lambda(p = .x$P.value))) %>%
  # Assign the model with inflation factor closes to 1 as the best model
  group_by(trait, parameter) %>%
  mutate(vif_best_model = nPC[which.min(abs(1 - p_value_infl))],
         vif_best_model = ifelse(nPC == vif_best_model, "*", "")) %>%
  ungroup()


# Choose the best model
univariate_gwas_farm_bestModel_vif <- univariate_gwas_farm_vif %>%
  filter(vif_best_model == "*") %>%
  select(trait, parameter, nPC)


## Examine QQ plots and determine the best model
# First write a CSV and then fill it in
univariate_gwas_farm_bestModel_vif %>%
  rename(bestModel_vif_nPC = nPC) %>%
  write_csv(x = ., file = file.path(result_dir, "univariate_gwas_farm_bestModel.csv"))
  
# Read in the filled version
univariate_gwas_farm_bestModel_qq <- read_csv(file = file.path(result_dir, "univariate_gwas_farm_bestModel_filled.csv"))




# Calculate variance explained --------------------------------------------

# M matrix
M <- t(as.matrix(GD)) - 1

univariate_gwas_farm_varexp <- univariate_gwas_farm %>%
  crossing(., kinship = c("farmcpu", "all_markers")) %>%
  rowid_to_column() %>%
  group_by(trait, parameter, nPC, kinship) %>%
  do({
    row <- .
    # print(row$rowid)
    kin_markers <- if (row$kinship == "farmcpu") row$kinship.markers[[1]] else colnames(M)
    response <- names(row$pred[[1]])[2]
    n_PC <- row$nPC
    
    # Calculate kinship matrix from the kin_markers
    if (!all(is.na(kin_markers))) {
      Km <- A.mat(X = M[levels(row$pred[[1]]$line_name), kin_markers, drop = TRUE], min.MAF = 0, max.missing = 1)
    } else {
      Km <- diag(ncol(Km))
    }
    
    
    # Significant markers
    gwas_sigmar <- as.character(subset(row$scores[[1]], P.value <= row$fdr_p_thresh, marker, drop = TRUE))
    
    # Combine PCs, sigmar, and Kinship
    mf <- cbind(row$pred[[1]], Kpca[levels(row$pred[[1]]$line_name), seq_len(n_PC), drop = FALSE], M[levels(row$pred[[1]]$line_name), gwas_sigmar, drop = FALSE])
    Z <- model.matrix(~ -1 + line_name, mf)
    X <- as.matrix(cbind(`(Intercept)` = 1, mf[,c(grep(pattern = "PC[0-9]{3}", x = names(mf), value = TRUE), gwas_sigmar), drop = FALSE]))
    
    # Fit a mixed model
    fit <- mixed.solve(y = mf[[response]], Z = Z, K = Km, X = X, method = "REML")
    
    # R2 of PCs, kinship, and significant markers
    y_hat_pc <- X[,grepl(pattern = "PC[0-9]{3}", x = colnames(X)), drop = FALSE] %*% fit$beta[grepl(pattern = "PC[0-9]{3}", x = names(fit$beta))]
    r2_pc <- cor(y_hat_pc, mf[[response]])^2
    r2_kin <- fit$Vu / var(mf[[response]])
    y_hat_sigmar <- X[, gwas_sigmar, drop = FALSE] %*% fit$beta[gwas_sigmar]
    r2_sigmar <- cor(y_hat_sigmar, mf[[response]])^2
    # R2 of each sigmar
    y_hat_sigmar_mat <- X[, gwas_sigmar, drop = FALSE] * matrix(fit$beta[gwas_sigmar], nrow = nrow(X), ncol = length(gwas_sigmar), byrow = TRUE)
    r2_sigmar_mat <- apply(X = y_hat_sigmar_mat, MARGIN = 2, FUN = function(x) cor(x, mf[[response]])^2)

    tibble(r2_PC = as.numeric(r2_pc), r2_kinship = as.numeric(r2_kin), r2_sigmar = as.numeric(r2_sigmar),
           r2_sigmar_df = list(tibble(marker = names(r2_sigmar_mat), r2 = as.numeric(r2_sigmar_mat))))
    
  }) %>% ungroup()








# Save --------------------------------------------------------------------




# Save everything
save("univariate_gwas_farm",  "univariate_gwas_farm_bestModel_qq", "univariate_gwas_farm_varexp",
     file = file.path(result_dir, "stability_gwas_out.RData"))





# 



