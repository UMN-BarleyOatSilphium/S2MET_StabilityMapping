## S2MET GWAS of Stability Coefficients
## 
## Perform association analysis of the FW stability coefficients for regular phenotype
## FW, and for FW on environmental covariables
## 
## Author: Jeff Neyhart
## Last updated: May 25, 2018
## 

library(qvalue)

# Run the source script
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))


# Extract the SNP information for later use
snp_info <- S2TP_imputed_multi_genos_hmp %>%
  select(marker = rs, chrom, pos, cM_pos)


# Load the FW results
load(file.path(result_dir, "S2MET_pheno_mean_fw_results.RData"))
# Load the FW sampling results
load(file.path(result_dir, "S2MET_pheno_fw_resampling.RData"))


## Association analysis
# Significance threshold
alpha <- 0.05

# Format the genotype data for use
geno_use <- S2TP_imputed_multi_genos_hmp %>%
  select(-alleles, -cM_pos) %>%
  rename(marker = rs) %>%
  as.data.frame()


# Format for modeling
pheno_to_model <- S2_MET_pheno_mean_fw %>% 
  distinct(line_name, trait, g, b, delta) %>% 
  mutate(log_delta = log(delta)) %>%
  select(-delta) %>%
  split(.$trait) %>%
  map(select, -trait) %>%
  map(as.data.frame)

# K matrix
M <- S2TP_imputed_multi_genos_mat
K <- A.mat(X = M, min.MAF = 0, max.missing = 1)


## Use the GWAS function in rrBLUP
## K and QK models
gwas_scan_QK <- pheno_to_model %>%
  map(~GWAS(pheno = ., geno = geno_use, K = K, n.PC = 1, min.MAF = 0, plot = FALSE))
gwas_scan_QK1 <- gwas_scan_QK %>% 
  list(., names(.)) %>% 
  pmap_df(~mutate(.x, trait = .y)) %>%
  mutate(model = "QK")

gwas_scan_K <- pheno_to_model %>%
  map(~GWAS(pheno = ., geno = geno_use, K = K, n.PC = 0, min.MAF = 0, plot = FALSE))
gwas_scan_K1 <- gwas_scan_K %>% 
  list(., names(.)) %>% 
  pmap_df(~mutate(.x, trait = .y)) %>%
  mutate(model = "K")

# Separate markers per chromosome
markers_per_chrom <- geno_use %>% split(.$chrom) %>% map("marker")
K_per_chrom <- markers_per_chrom %>% map(~setdiff(geno_use$marker, .)) %>% 
  map(~A.mat(X = S2TP_imputed_multi_genos_mat[,.], min.MAF = 0, max.missing = 1))

# QG model
gwas_scan_QG <- pheno_to_model %>%
  map(function(p) {
    
    geno_use %>% 
      split(.$chrom) %>%
      list(., K_per_chrom) %>%
      pmap(~GWAS(pheno = p, geno = .x, K = .y, n.PC = 1, min.MAF = 0, plot = FALSE)) %>%
      bind_rows()
    
  })

gwas_scan_QG1 <- gwas_scan_QG %>% 
  list(., names(.)) %>% 
  pmap_df(~mutate(.x, trait = .y)) %>%
  mutate(model = "QG")

# QG model
gwas_scan_G <- pheno_to_model %>%
  map(function(p) {
    
    geno_use %>% 
      split(.$chrom) %>%
      list(., K_per_chrom) %>%
      pmap(~GWAS(pheno = p, geno = .x, K = .y, n.PC = 0, min.MAF = 0, plot = FALSE)) %>%
      bind_rows()
    
  })

gwas_scan_G1 <- gwas_scan_G %>% 
  list(., names(.)) %>% 
  pmap_df(~mutate(.x, trait = .y)) %>%
  mutate(model = "G")

## Combine all the results
gwas_pheno_mean_fw <- bind_rows(gwas_scan_K1, gwas_scan_G1, gwas_scan_QK1, gwas_scan_QG1) %>%
  select(trait, model, names(.))

# Tidy
gwas_pheno_mean_fw_tidy <- gwas_pheno_mean_fw %>% 
  gather(coef, neg_log_p, g, b, log_delta)


## Save the data
save("gwas_pheno_mean_fw_tidy", file = file.path(result_dir, "pheno_fw_mean_gwas_results.RData"))





## Load data
load(file.path(result_dir, "pheno_fw_mean_gwas_results.RData"))


## Plot the QQ plot
gwas_pheno_mean_fw_qq <- gwas_pheno_mean_fw_tidy %>% 
  group_by(trait, model, coef) %>% 
  arrange(trait, model, coef, desc(neg_log_p)) %>% 
  mutate(neg_log_p_exp = -log10(ppoints(n = n())),
         # Add a confidence interval based on the beta distribution (assumes independence of tests)
         ci_lower = -log10(qbeta(p = (alpha / 2), shape1 = seq(n()), rev(seq(n())))),
         ci_upper = -log10(qbeta(p = 1 - (alpha / 2), shape1 = seq(n()), rev(seq(n())))))

g_gwas_qq <- gwas_pheno_mean_fw_qq %>% 
  ggplot(aes(x = neg_log_p_exp, y = neg_log_p, color = model)) + 
  geom_ribbon(aes(ymin = ci_upper, ymax = ci_lower), fill = "grey75") +
  geom_point() + 
  geom_abline(aes(slope = 1, intercept = 0)) +
  facet_grid(coef ~ trait) +
  theme_bw()

# Save
save_file <- file.path(fig_dir, "pheno_fw_mean_gwas_qq.jpg")
ggsave(filename = save_file, plot = g_gwas_qq, width = 8, height = 8, dpi = 1000)


# Convert p-value to q-values
gwas_pheno_mean_fw_tidy_adj <- gwas_pheno_mean_fw_tidy %>% 
  group_by(trait, model, coef) %>% 
  mutate(qvalue = qvalue(p = 10^-neg_log_p)$qvalue, 
         neg_log_q = -log10(qvalue),
         color = if_else(chrom %in% seq(1, 7, 2), "B", "G")) %>%
  ungroup()
           

# Plot modifier
# Manhattan plot
g_mod_man <- list(
  geom_point(),
  geom_hline(yintercept = -log10(alpha), lty = 2),
  # geom_hline(aes(yintercept = neg_log10_fdr10, lty = "FDR 10%")),
  scale_color_manual(values = color, guide = FALSE),
  ylab(expression(-log[10](italic(q)))),
  xlab("Position (Mbp)"),
  theme_bw(),
  theme_manhattan(),
  theme(panel.border = element_blank())
)

# Iterate over the models and plot
for (mod in unique(gwas_pheno_mean_fw_tidy_adj$model)) {
  
  # Create the plot
  g_gwas_manhattan <- gwas_pheno_mean_fw_tidy_adj %>%
    filter(model == mod) %>%
    ggplot(aes(x = pos / 1e6, y = neg_log_q, group = chrom, col = color)) + 
    facet_grid(trait + coef ~ chrom, switch = "x", scales = "free", space = "free_x") +
    g_mod_man
  
  ggsave(filename = str_c("gwas_manhattan_allcoef_", mod, ".jpg"), plot = g_gwas_manhattan,
         path = fig_dir, height = 10, width = 8, dpi = 1000)
  
}


## Iterate over traits and plot
for (tr in unique(gwas_pheno_mean_fw_tidy_adj$trait)) {
  
  # Create the plot
  g_gwas_manhattan <- gwas_pheno_mean_fw_tidy_adj %>%
    filter(trait == tr) %>%
    ggplot(aes(x = pos / 1e6, y = neg_log_q, group = chrom, col = color)) + 
    facet_grid(coef + model ~ chrom, switch = "x", scales = "free", space = "free_x") +
    g_mod_man
  
  ggsave(filename = str_c("gwas_manhattan_allcoef_allmod_", tr, ".jpg"), plot = g_gwas_manhattan,
         path = fig_dir, height = 12, width = 8, dpi = 1000)
  
}







## Fit a multi-locus mixed model to identify significant QTL

# Identify the signficant associations
gwas_adj_sig <- gwas_pheno_mean_fw_tidy_adj %>%
  filter(qvalue <= alpha)



# Adjust the phenotypic data
pheno_to_model_mlmm <- pheno_to_model %>%
  list(., names(.)) %>%
  pmap_df(~mutate(.x, trait = .y)) %>%
  gather(coef, value, g:log_delta)



# Iterate over trait, coef, and chromosome
gwas_mlmm_model <- gwas_adj_sig %>%
  group_by(trait, coef, model) %>%
  do(mlmm_out = {
    # Extract the data
    df <- .
    
    # What are the traits and coefficients?
    tr <- unique(df$trait)
    cf <- unique(df$coef)
    mod <- unique(df$model)

    # cat(tr, "\t", cf, "\t", mod, "\n")
    
    # Subset the data.frame and create a model.frame
    mf <- pheno_to_model_mlmm %>%
      filter(trait == tr, coef == cf) %>%
      model.frame(formula = value ~ line_name)

    # Response vector
    y <- model.response(mf)
    # Fixed effect mean matrix
    X_mu <- model.matrix(~ 1, mf)
    # Random effect of each entry
    Z_g <- model.matrix(~ -1 + line_name, mf)
    # Fixed effect of each SNP
    X_snp <- remove_colinear((Z_g %*% M[,df$marker, drop = FALSE]))
    
    # If the model contains Q, calculate pop structure
    if (grepl("Q", mod)) {
      Q <- eigen(x = K)$vector[,1]
    } else {
      Q <- NULL
    }

    # Combine the fixed effect matrices
    X <- cbind(X_mu, X_snp, Q = Q)

    # Fit the model and return estimates
    fit_out <- fit_out_base <- mlmm(y = y, X = X, Z = Z_g, K = K, snps = colnames(X_snp))

    # Are any qvalues greater than the significance threshold?
    snp_summary <- subset(fit_out$summary, term %in% colnames(X_snp)) 
    nonsig <- snp_summary$qvalue > alpha

    # While loop to backwards eliminate
    while(any(nonsig)) {

      # Which markers had the largest p_value
      which_remove <- which.max(snp_summary$qvalue)
      # Remove that marker from the X matrix
      X_snp <- X_snp[,-which_remove, drop = FALSE]

      # If no SNPs are present, return NA
      if (ncol(X_snp) == 0) {
        fit_out <- list(NA)
        # Stop the loop
        break

      } else {
        # New X matrix
        X_new <- cbind(X_mu, X_snp, Q = Q)

        # Fit the model and return estimates
        fit_out <- mlmm(y = y, X = X_new, Z = Z_g, K = K, snps = colnames(X_snp))

      }

      # Are any qvalues greater than the significance threshold?
      snp_summary <- subset(fit_out$summary, term %in% colnames(X_snp)) 
      nonsig <- snp_summary$qvalue > alpha

    } # End of while loop
    
    # Return a list
    list(fit_out_base = fit_out_base, fit_out_reduced = fit_out)
    
  }) %>% ungroup()


## Get the estimates of SNP effects
gwas_mlmm_final <- gwas_mlmm_model %>% 
  mutate(mlmm_reduced = map(mlmm_out, "fit_out_reduced"),
         mlmm_reduced = map(mlmm_reduced, ~cbind(subset(.$summary, term != "Q"), r_squared = .$r_squared$snp_r_squared, 
                                                 t(.$r_squared$fixed_r_squared))) ) %>%
  unnest(mlmm_reduced)


# Save this data
save_file <- file.path(result_dir, "pheno_fw_mean_gwas_results.RData")
save("gwas_pheno_mean_fw_tidy_adj", "gwas_adj_sig", "gwas_mlmm_model",
     "gwas_mlmm_model", file = save_file)





### Mapping using resampling estimates of stability ###


## Use the FW resampling data to determine the robustness of the mapping results
# Tidy up for splitting by cores
resample_phenos_use <- S2MET_pheno_sample_fw %>%
  # Convert delta to log_delta
  mutate(log_delta = log(delta)) %>% 
  # Tidy
  select(-delta) %>% 
  group_by(trait, p, iter) %>% 
  nest()

# Assign cores
resample_phenos_use_list <- resample_phenos_use %>%
  ungroup() %>%
  mutate(core = sort(rep(seq(n_core), length.out = nrow(.)))) %>%
  split(.$core)


## Parallelize the association
resample_gwas_sig_out <- mclapply(X = resample_phenos_use_list, function(core_df) {
  
  # Map over the list of data
  gwas_out_list <- core_df$data %>%
    map(~gwas(pheno = ., geno = geno_use, model = "G", P3D = TRUE, n.core = 1))
  
  # Grab the scores, adjust the pvalues, then filter for those <= alpha
  gwas_sig_scores <- gwas_out_list %>% 
    map("scores") %>%
    map(~mutate(., padj = p.adjust(pvalue, method = "fdr"))) %>%
    map(~filter(., padj <= alpha))
  
  # Add those scores to the data.frame and return
  core_df %>% 
    select(-data, -core) %>% 
    mutate(sig_scores = gwas_sig_scores)
  
}, mc.cores = n_core)





# Save the results
save_file <- file.path(result_dir, str_c("S2MET_pheno_fw_gwas_resample_results_", tr, ".RData"))
save("resample_gwas_sig_out", file = save_file)





