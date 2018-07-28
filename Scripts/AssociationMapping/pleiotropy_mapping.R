## S2MET GWAS for pleiotropy between genotype mean and stability
## 
## Perform pleiotropic association analysis of the FW stability coefficients and 
## genotype means
## 
## Author: Jeff Neyhart
## Last updated: May 30, 2018
## 


# Run the source script - local
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

library(EMMREML)
library(qvalue)
library(parallel)


# Load the FW results
load(file.path(result_dir, "pheno_mean_fw_results.RData"))
# Load the FW sampling results
load(file.path(result_dir, "pheno_fw_resampling.RData"))

## Number of cores
n_core <- ifelse(Sys.info()["sysname"] == "Windows", 1, detectCores())


# Format the genotype data for use
geno_use <- s2tp_genos_imputed_hmp %>%
  select(-alleles, -cM_pos) %>%
  as.data.frame()

pheno_mean_fw <- pheno_mean_fw_tpvp %>%
  filter(line_name %in% tp)


# Format for modeling
pheno_to_model <- pheno_mean_fw %>% 
  distinct(line_name, trait, g, b, delta) %>% 
  mutate(log_delta = log(delta)) %>%
  select(-delta) %>%
  split(.$trait) %>%
  map(select, -trait) %>%
  map(as.data.frame)

# K matrix
M <- s2tp_genos_imputed
K <- A.mat(X = M, min.MAF = 0, max.missing = 1)
# Q vector
Q <- eigen(x = K)$vector[,1]

# K matrix by chromosome
markers_per_chrom <- geno_use %>% 
  split(.$chrom) %>% 
  map("marker")
K_per_chrom <- markers_per_chrom %>% 
  map(~setdiff(geno_use$marker, .)) %>% 
  map(~A.mat(X = s2tp_genos_imputed[,.], min.MAF = 0, max.missing = 1))
Q_per_chrom <- map(K_per_chrom, ~eigen(x = .)$vector[,1])



## Models to fit
# models <- c("K", "QK", "G", "QG")
models <- "QK"


## Test for pleiotropy
##
## # Fit a multi-variate mixed model for the genotype mean and stability,
## # then test each marker using an intersection-union test
##
## # H0: the marker is not associated with the mean or it is not associated with stability
## # HA: the marker is associated with both mean and stability
##

# Re-organize the data for modeling
pheno_to_model_plei <- pheno_to_model %>% 
  list(., names(.)) %>% 
  pmap_df(~mutate(., trait = .y)) %>% 
  gather(stab_coef, value, b, log_delta)

# Split up the data for iterations
pheno_to_model_plei_split <- pheno_to_model_plei %>%
  filter(stab_coef == "b") %>%
  split(list(.$trait, .$stab_coef)) %>%
  map(~spread(., stab_coef, value) %>% select(-trait))

# ## Run the GWAS for pleiotropy
# # K and QK models
# gwas_mv_scan_QK <- pheno_to_model_plei_split %>%
#   map(~mv_gwas(pheno = ., geno = geno_use, K = K, n.PC = 1))
# gwas_mv_scan_QK1 <- gwas_mv_scan_QK %>% 
#   list(., names(.)) %>% 
#   pmap_df(~mutate(.x, trait = .y)) %>%
#   mutate(model = "QK", trait = str_replace(string = trait, pattern = "\\.log_delta|\\.b", ""))

# gwas_mv_scan_K <- pheno_to_model_plei_split %>%
#   map(~mv_gwas(pheno = ., geno = geno_use, K = K, n.PC = 0))
# gwas_mv_scan_K1 <- gwas_mv_scan_K %>% 
#   list(., names(.)) %>% 
#   pmap_df(~mutate(.x, trait = .y)) %>%
#   mutate(model = "K", trait = str_replace(string = trait, pattern = "\\.log_delta|\\.b", ""))


# QG and G models
gwas_mv_scan_QG <- pheno_to_model_plei_split %>%
  map(function(p) {
    
    geno_use %>% 
      split(.$chrom) %>%
      list(., K_per_chrom) %>%
      pmap(~mv_gwas(pheno = p, geno = .x, K = .y, n.PC = 1)) %>%
      bind_rows()
    
  })

gwas_mv_scan_QG1 <- gwas_mv_scan_QG %>% 
  list(., names(.)) %>% 
  pmap_df(~mutate(.x, trait = .y)) %>%
  mutate(model = "QG", trait = str_replace(string = trait, pattern = "\\.log_delta|\\.b", ""))


# gwas_mv_scan_G <- pheno_to_model_plei_split %>%
#   map(function(p) {
#     
#     geno_use %>% 
#       split(.$chrom) %>%
#       list(., K_per_chrom) %>%
#       pmap(~mv_gwas(pheno = p, geno = .x, K = .y, n.PC = 0)) %>%
#       bind_rows()
#     
#   })
# 
# gwas_mv_scan_G1 <- gwas_mv_scan_G %>% 
#   list(., names(.)) %>% 
#   pmap_df(~mutate(.x, trait = .y)) %>%
#   mutate(model = "G", trait = str_replace(string = trait, pattern = "\\.log_delta|\\.b", ""))

## Combine all the results
gwas_pleio_pheno_mean_fw <- bind_rows(
  # gwas_mv_scan_K1, 
  # gwas_mv_scan_G1, 
  # gwas_mv_scan_QK1, 
  gwas_mv_scan_QG1) %>%
  select(trait, model, names(.))



## Save the data
save_file <- file.path(result_dir, "pheno_fw_gwas_pleiotropy_results.RData")
save("gwas_pleio_pheno_mean_fw", file = save_file)




## Load the data
load(file.path(result_dir, "pheno_fw_gwas_pleiotropy_results.RData"))

# ## Run a shortcut of the intersection-union test by grabbing the minimum -log10(p)
# gwas_pleio_pheno_mean_fw1 <- gwas_pleio_pheno_mean_fw %>% 
#   mutate(coef_pair = pmap(list(g, b, log_delta), ~which(!is.na(c(..1, ..2, ..3)))) %>% 
#            map(~c("g", "b", "log_delta")[.]) %>% map_chr(~paste0(., collapse = ":")),
#          min_neg_log_p = pmap_dbl(list(g, b, log_delta), min, na.rm = T)) %>%
#   select(-g:-log_delta)


gwas_pleio_pheno_mean_fw1 <- gwas_pleio_pheno_mean_fw %>%
  gather(coef_pair, min_neg_log_p, -trait:-pos)



## Plot the QQ plot
gwas_pleio_pheno_mean_fw_qq <- gwas_pleio_pheno_mean_fw1 %>% 
  group_by(trait, model, coef_pair) %>% 
  arrange(trait, model, coef_pair, desc(min_neg_log_p)) %>% 
  mutate(neg_log_p_exp = -log10(ppoints(n = n())),
         neg_log_p_exp_pair = -log10(sqrt(ppoints(n = n()))),
         # Add a confidence interval based on the beta distribution (assumes independence of tests)
         ci_lower = -log10(qbeta(p = (alpha / 2), shape1 = seq(n()), rev(seq(n())))),
         ci_upper = -log10(qbeta(p = 1 - (alpha / 2), shape1 = seq(n()), rev(seq(n())))),
         ci_lower_pair = -log10(sqrt(qbeta(p = (alpha / 2), shape1 = seq(n()), rev(seq(n()))))),
         ci_upper_pair = -log10(sqrt(qbeta(p = 1 - (alpha / 2), shape1 = seq(n()), rev(seq(n())))))) %>%
  ungroup()


g_gwas_pleio_qq <- gwas_pleio_pheno_mean_fw_qq %>% 
  # filter(trait == "GrainYield") %>%
  gather(test, expected_neg_log_p, neg_log_p_exp, neg_log_p_exp_pair) %>%
  mutate(test = ifelse(test == "neg_log_p_exp", "alpha", "sqrt(alpha)")) %>%
  ggplot(aes(x = expected_neg_log_p, y = min_neg_log_p, color = model, shape = test)) + 
  # geom_ribbon(aes(ymin = ci_upper, ymax = ci_lower), fill = "grey75") +
  # geom_ribbon(aes(ymin = ci_upper_pair, ymax = ci_lower_pair), fill = "grey75") +
  geom_point() + 
  geom_abline(aes(slope = 1, intercept = 0)) +
  facet_grid(coef_pair ~ trait) +
  theme_bw()





# Convert p-value to q-values
gwas_pleio_pheno_mean_fw_tidy_adj <- gwas_pleio_pheno_mean_fw1 %>% 
  group_by(trait, model, coef_pair) %>% 
  mutate(qvalue = qvalue(p = 10^-min_neg_log_p)$qvalue, 
         neg_log_q = -log10(qvalue),
         color = if_else(chrom %in% seq(1, 7, 2), "B", "G"),
         empirical_cutoff = quantile(x = 10^-min_neg_log_p, probs = alpha / 2)) %>%
  ungroup()



# Plot modifier
# Manhattan plot
g_mod_man <- list(
  geom_point(),
  # geom_hline(yintercept = -log10(alpha), lty = 2),
  geom_hline(aes(yintercept = -log10(empirical_cutoff)), lty = 2),
  # geom_hline(aes(yintercept = neg_log10_fdr10, lty = "FDR 10%")),
  scale_color_manual(values = color, guide = FALSE),
  ylab(expression(-log[10](italic(p)))),
  xlab("Position (Mbp)"),
  theme_bw(),
  theme_manhattan(),
  theme(panel.border = element_blank())
)

# Iterate over the models and plot
for (mod in unique(gwas_pleio_pheno_mean_fw_tidy_adj$model)) {
  
  # Create the plot
  g_gwas_manhattan <- gwas_pleio_pheno_mean_fw_tidy_adj %>%
    filter(model == mod) %>%
    ggplot(aes(x = pos / 1e6, y = min_neg_log_p, group = chrom, col = color)) +
    # ggplot(aes(x = pos / 1e6, y = neg_log_q, group = chrom, col = color)) + 
    facet_grid(trait + coef_pair ~ chrom, switch = "x", scales = "free", space = "free_x") +
    g_mod_man
  
  ggsave(filename = str_c("gwas_pleio_manhattan_allcoef_", mod, ".jpg"), plot = g_gwas_manhattan,
         path = fig_dir, height = 10, width = 8, dpi = 1000)
  
}


