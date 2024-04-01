## Analysis of phenotypic stability
## 
## Ranges, plots, variance components
## 
## Author: Jeff Neyhart
## Last updated: June 23, 2018
## 

# Other packages
library(cowplot)
library(modelr)
library(EMMREML)
library(lme4qtl)

# Repository directory
repo_dir <- getwd()
# Project and other directories
source(file.path(repo_dir, "source.R"))
 
# Load fw data
load(file.path(result_dir, "pheno_mean_fw_results.RData"))

# Relationship matrix
M <- s2tp_genos_imputed
K <- A.mat(X = M, min.MAF = 0, max.missing = 1)

# Significance threshold
alpha <- 0.05

# # Substitute the estimates of stability using just the TP with estimates using
# # both the TP and VP
# pheno_mean_fw <- pheno_mean_fw_tpvp %>%
#   filter(line_name %in% tp)

pheno_mean_fw1 <- pheno_mean_fw %>%
  filter(set == "herit_environments") %>%
  bind_rows(., filter(pheno_mean_fw, set != "herit_environments", !trait %in% unique(.$trait)))


### Preliminary analysis

## Are traits correlated with latitude?
env_means <- pheno_mean_fw1 %>% 
  distinct(set, trait, environment, h) %>%
  left_join(distinct(trial_info, environment, latitude, longitude))

env_means %>%
  group_by(set, trait) %>%
  summarize_at(vars(latitude, longitude), funs(list(cor.test(., h)))) %>%
  gather(coord, out, -trait, -set) %>%
  mutate(cor = map_dbl(out, "estimate"),
         p_value = map_dbl(out, "p.value")) %>%
  filter(set == "herit_environments")





### Finlay-Wilkinson Regression

## Results

# What are the ranges of genotype values and environmental effects per trait?
pheno_summary <- pheno_mean_fw1 %>%
  distinct(trait, line_name, environment, h, g)  %>% 
  gather(term, value, g:h) %>% 
  group_by(trait, term) %>% 
  summarize_at(vars(value), funs(min, max)) %>% 
  mutate(range = max - min)

## Environment metadata
trial_info %>% 
  filter(environment %in% pheno_mean_fw1$environment) %>%
  select(environment, latitude, longitude) %>%
  gather(measure, value, -environment) %>% 
  group_by(measure) %>% 
  summarize_at(vars(value), funs(min, max), na.rm = T) %>% 
  mutate(range = max - min)


### Phenotypic Stability

# Transformation

# Extract the unique stability coefficients
# Then add the breeding program information
pheno_mean_fw2 <- pheno_mean_fw1 %>%
  left_join(., subset(entry_list, Class == "S2TP", c(Line, Program)), 
            by = c("line_name" = "Line")) %>%
  rename(program = Program)

## Log transform the non-linear stability estimates
pheno_fw_use <- pheno_mean_fw2 %>%
  group_by(trait) %>% 
  mutate(log_delta = log(delta)) %>%
  # Tidy
  gather(term, estimate, b, delta, log_delta) %>% 
  filter(term != "delta")


#### Summary

# Ranges and variances

# What is the range of linear stability per trait?
pheno_fw_use %>%
  filter(term == "b") %>%
  group_by(trait) %>% 
  summarize_at(vars(estimate), funs(min, max)) %>% 
  mutate(range = max - min)

# Extract the distinct combinations of line, trait, and stability
pheno_fw_uniq <- pheno_fw_use %>%
  distinct(trait, line_name, g, program, term, estimate)


# Plot the typical FW regression (i.e. using the phenotypic means)

## Plot the stability terms individually, then combine
# First define a common list of ggplot modifiers
g_mod <- list(
  geom_density(aes(fill = "blue")),
  xlab("Estimate"),
  theme_pnas() +
  theme(axis.title.y = element_blank()) )

# Just plot linear stability
g_pheno_fw_dens_b <- pheno_fw_uniq %>%
  filter(term == "b") %>%
  ggplot(aes(x = estimate)) + 
  labs(title = expression("Linear Stability"~(b[i]))) +
  facet_wrap( ~ trait, ncol = 1) +
  scale_fill_brewer(palette = "Set2", guide = FALSE) +
  g_mod


# Just plot non-linear stability
g_pheno_fw_dens_delta <- pheno_fw_uniq %>%
  filter(term == "log_delta") %>%
  ggplot(aes(x = estimate)) + 
  labs(title = expression("Non-Linear Stability"~(ln(delta[i]^2)))) +
  scale_fill_brewer(palette = "Set2", guide = FALSE) +
  facet_wrap( ~ trait, ncol = 1, scales = "free_x") +
  g_mod

# Add the plots together
g_fw_dist <- plot_grid(g_pheno_fw_dens_b, g_pheno_fw_dens_delta, ncol = 2, align = "hv")

ggsave(filename = "stability_estimate_distriubtions.jpg", plot = g_fw_dist, path = fig_dir,
       height = 10, width = 5, dpi = 1000)


# Just plot linear stability
g_pheno_fw_dens_b <- pheno_fw_uniq %>%
  filter(term == "b", trait %in% traits_specific) %>%
  ggplot(aes(x = estimate)) + 
  labs(title = expression("Linear Stability"~(b[i]))) +
  facet_wrap( ~ trait, ncol = 1) +
  scale_fill_brewer(palette = "Set2", guide = FALSE) +
  g_mod


# Just plot non-linear stability
g_pheno_fw_dens_delta <- pheno_fw_uniq %>%
  filter(term == "log_delta", trait %in% traits_specific) %>%
  ggplot(aes(x = estimate)) + 
  labs(title = expression("Non-Linear Stability"~(ln(delta[i]^2)))) +
  scale_fill_brewer(palette = "Set2", guide = FALSE) +
  facet_wrap( ~ trait, ncol = 1, scales = "free_x") +
  g_mod

# Add the plots together
g_fw_dist <- plot_grid(g_pheno_fw_dens_b, g_pheno_fw_dens_delta, ncol = 2, align = "hv")


ggsave(filename = "stability_estimate_distriubtions_specific.jpg", plot = g_fw_dist, path = fig_dir,
       height = 10, width = 8.7, units = "cm", dpi = 1000)




## What is the proportion of variation in the mean and stability expained by genomewide markers?
pheno_fw_use_tomodel <- pheno_fw_use %>%
  filter(trait %in% traits_specific) %>%
  distinct(trait, line_name, g, term, estimate) %>%
  spread(term, estimate) %>%
  gather(term, estimate, g:log_delta) %>%
  filter(line_name %in% tp_geno) %>%
  mutate(markers = line_name)

# Group by trait and term
pheno_fw_marker_varcomp <- pheno_fw_use_tomodel %>%
  group_by(trait, term) %>%
  do({
    df <- .
    # Fit a model
    fit <- relmatLmer(estimate ~ (1|markers) + (1|line_name), data = df, relmat = list(markers = K))
    
    # Return the variance components and proportions
    VarProp(fit) %>% 
      select(grp, vcov) %>% 
      add_row(grp = "total_genetic", vcov = sum(.$vcov[1:2])) %>% 
      mutate(prop_total = ifelse(grp == "total_genetic", NA, vcov / sum(head(vcov, -1))), 
             prop_genetic = ifelse(grp == "Residual", NA, vcov / tail(vcov, 1)))
    
  })

## Output a table
write_csv(x = pheno_fw_marker_varcomp, path = file.path(fig_dir, "trait_marker_varcomp.csv"))















# Plot the stability estimates as lines against g_i vs h_j
# colors_use <- set_names(umn_palette(2)[3:7], unique(pheno_fw_use$program))
color_use <- umn_palette(2)[3]

# Add the program information to the results
pheno_fw_use_toplot <- pheno_fw_use %>% 
  filter(trait %in% traits_specific) %>%
  spread(term, estimate) %>% 
  select(trait, environment, line_name, value, g, h, b, log_delta) %>% 
  left_join(., subset(entry_list, Class == "S2TP", c(Line, Program)), by = c("line_name" = "Line")) %>%
  rename(program = Program) %>%
  ungroup()

# Plot the normal FW analysis plot (genotype mean against environmental effects)
g_pheno_fw_b <- pheno_fw_use_toplot %>%
  ggplot(aes(x = h, y = value, group = line_name)) + 
  geom_point(size = 0.1, color = "black") + 
  geom_abline(aes(slope = b, intercept = g), color = color_use, alpha = 0.15) +
  facet_wrap(~ trait, ncol = 1, scales = "free", strip.position = "left") +
  # scale_color_manual(drop = FALSE, name = "Origin\nBreeding\nProgram", values = colors_use,
  #                    guide = guide_legend(override.aes = list(size = 1), nrow = 1)) +
  ylab("Phenotypic value") +
  xlab("Environmental mean") +
  # labs(title = "Linear Phenotypic Stability") +
  # scale_color_brewer(palette = "Set2") +
  theme_pnas() +
  theme(legend.position = "bottom", legend.direction = "horizontal", strip.placement = "outside")

# ggsave(filename = "pheno_fw_linear.jpg", plot = g_pheno_fw_b, path = fig_dir, 
#        width = 5, height = 10, dpi = 1000)


# ## Subset genotypes to highlight different responses
# pheno_fw_example <- pheno_fw_use_toplot %>%
#   distinct(trait, line_name, g, b, program) %>%
#   mutate(program = as.factor(program)) %>%
#   group_by(trait) %>% 
#   filter(b == min(b) | b == max(b) | abs(1 - b) == min(abs(1 - b))) %>%
#   left_join(., data.frame(program = unique(pheno_fw_use_toplot$program), color = colors_use))


# g_pheno_fw_b_example <- pheno_fw_use_toplot %>%
#   select(trait, environment, line_name, value, h) %>% 
#   ggplot(aes(x = h, y = value)) + 
#   geom_point(color = "grey75", size = 0.2) + 
#   geom_abline(data = pheno_fw_example, aes(intercept = g, slope = b, col = program), lwd = 0.5) + 
#   # scale_color_manual(drop = FALSE, guide = guide_legend(title = "Program", ncol = 1), values = colors_use) +
#   scale_color_manual(drop = FALSE, guide = FALSE, values = colors_use) +
#   facet_wrap(~ trait, ncol = 1, scales = "free", strip.position = "left") +
#   ylab("Phenotypic Value") +
#   xlab("Environment Mean") +
#   # labs(title = "Example Genotype Responses") +
#   theme_pnas() +
#   theme(legend.position = c(0.10, 0.90), axis.title.y = element_blank(),
#         strip.background = element_blank(), strip.text = element_blank())
# 
# ggsave(filename = "pheno_fw_linear_example.jpg", plot = g_pheno_fw_b_example, path = fig_dir, 
#        width = 5, height = 10, dpi = 1000)







## Plot the correlation between the genotypic effect and the sensitivity
set.seed(415)

n_perm <- 5000

## Use permutation to estimate the correlation
stability_mean_corr_perm <- pheno_fw_use %>%
  distinct(trait, line_name, g, term, estimate) %>%
  filter(!is.infinite(estimate), trait %in% traits_specific) %>%
  group_by(trait, term) %>%
  do({
    df <- .
    
    # Estimate the base correlation
    base_cor <- cor(df$g, df$estimate)
    
    perms <- permute(data = df, n = n_perm, estimate)
    perm_cor <- map_dbl(perms$perm, ~{
      d <- as.data.frame(.)
      cor(d$g, d$estimate)
    })
    
    # Estimate the p-value
    pvalue <- mean(perm_cor >= abs(base_cor) | perm_cor <= -abs(base_cor))
    
    # Confidence intervals
    ci <- quantile(perm_cor, c(alpha / 2, 1 - (alpha / 2))) + base_cor
    ci_lower <- ci[1]
    ci_upper <- ci[2]
    
    data_frame(base = base_cor, pvalue = pvalue, ci_lower = ci_lower, ci_upper = ci_upper)
    
  }) %>%
  # Add annotation
  mutate(significant = !between(0, ci_lower, ci_upper),
         sig_ann = ifelse(significant, "*", ""),
         annotation = str_c("r = ", round(base, 3), sig_ann))
  


# Create a list of plot additions
g_add <- list(geom_point(size = 0.1, color = "black"),
              geom_smooth(method = "lm", se = FALSE, col = color_use, lwd = 0.5),
              # scale_color_manual(name = "Program", values = colors_use, guide = FALSE),
              xlab("Genotype mean"),
              theme_pnas(),
              theme(strip.placement = "inside"))

# Plot just the linear stability
g_linear_stability_and_mean <- pheno_fw_use %>%
  filter(term == "b", trait %in% traits_specific) %>%
  ggplot(aes(x = g, y = estimate, group = FALSE)) +
  geom_text(data = subset(stability_mean_corr_perm, term == "b"), 
            aes(x = Inf, y = -Inf, label = annotation, vjust = -1, hjust = 1.2), col = "black", size = 1.5) +
  facet_wrap(~ trait, scales = "free_x", ncol = 1, strip.position = "left") +
  ylab("Linear stability estimate") +
  theme(legend.position = "right") +
  g_add

# Plot just the non-linear stability
g_nonlinear_stability_and_mean <- pheno_fw_use %>%
  filter(term == "log_delta", trait %in% traits_specific) %>%
  ggplot(aes(x = g, y = estimate, group = FALSE)) +
  geom_text(data = subset(stability_mean_corr_perm, term == "log_delta"),
            aes(x = Inf, y = -Inf, label = annotation, vjust = -1, hjust = 1.2), col = "black", size = 1.5) +
  facet_wrap(~ trait, scales = "free", ncol = 1, strip.position = "left")+
  ylab("Non-linear stability estimate") +
  theme(legend.position = "right") +
  g_add


# Create a plot of just the strip
g_strips <- g_pheno_fw_b + theme(legend.position = "none", strip.placement = "outside", axis.title = element_blank(), 
                                 axis.text = element_blank(), axis.ticks = element_blank(), axis.line = element_blank(), 
                                 panel.ontop = T)


## Plot the reaction norms, then the correlations together
g_plotgrid <- plot_grid(
  g_strips,
  g_pheno_fw_b + theme(legend.position = "none", strip.placement = "inside", strip.background = element_blank(), 
                       strip.text = element_blank()),
  g_linear_stability_and_mean + theme(strip.background = element_blank(), strip.text = element_blank()),
  g_nonlinear_stability_and_mean + theme(strip.background = element_blank(), strip.text = element_blank()),
  nrow = 1, align = "hv", labels = c("", LETTERS[1:3]), axis = "tblr", label_size = 8, rel_widths = c(0.15, 1, 1, 1)
)



# Add the legend
# g_plotgrid1 <- plot_grid(g_plotgrid, get_legend(g_pheno_fw_b), ncol = 1, rel_heights = c(1, 0.1))


ggsave(filename = "reaction_norms_and_correlations.jpg", plot = g_plotgrid, path = fig_dir, 
       width = 11.4, height = 8, units = "cm", dpi = 1000)





### Use a bi-variate model and the genomic relationship matrix to estimate genetic correlation
# Load permutation results
load(file.path(result_dir, "stability_correlation_permutation.RData"))

# Create a df for modeling
pheno_fw_use_tomodel <- pheno_fw_use %>%
  filter(trait %in% traits_specific,
         line_name %in% tp_geno) %>%
  ungroup() %>%
  distinct(trait, line_name, g, term, estimate)

# Iterate over each trait and term
pheno_fw_gen_corr1 <- pheno_fw_use_tomodel %>%
  group_by(trait, term) %>%
  do({
    df <- .
    # Create a model.frame
    mf <- model.frame(estimate ~ line_name, df)
    
    # Create model matrices
    Z <- model.matrix(~ -1 + line_name, mf) %>%
      `colnames<-`(., colnames(K))
    X <- model.matrix(~ 1, mf)
    
    # Create the response matrix
    Y <- as.matrix(select(df, g, estimate))
    # Y <- scale(Y)
    
    # Fit the model
    fit <- emmremlMultivariate(Y = t(Y), X = t(X), Z = t(Z), K = K)
    vcovG <- fit$Vg
    
    # # Use sommer
    # fit <- sommer::mmer(Y = Y, X = X, Z = list(gen = list(Z = Z, K = K)))
    # vcovG <- fit$var.comp$gen
    
    # Estimate the correlation using the formula for correlation
    # covariance / (sd1 * sd2)
    rhoG_hat <- vcovG[1,2] / prod(sqrt(diag(vcovG)))
    
    # Return
    data_frame(trait1 = "g", trait2 = unique(df$term), corG = rhoG_hat, df = nrow(df) - 2)
  })


# Two-tailed test using permutation results
pheno_fw_gen_corr_perm %>%
  filter(!is.na(corrG)) %>%
  left_join(., pheno_fw_gen_corr) %>% 
  group_by(trait, term) %>% 
  summarize(p_value = mean(corrG >= abs(corG) | corrG <= -abs(corG)),
            p_value_upper = mean(corrG >= corG),
            p_value_lower = mean(corrG <= corG))





## Calculate the genetic correlation after correcting for significant GWAS markers
# Load the GWAS results
load(file.path(result_dir, "pheno_fw_mean_gwas_results.RData"))


## Adjust the mlmm results
gwas_mlmm_model_adj <- subset(gwas_mlmm_final, model == "QG") %>%
  left_join(., select(snp_info, marker:pos), c("term" = "marker")) %>%
  select(set, trait:model, marker = term, chrom:pos, beta, pvalue, qvalue, snp_r_squared = r_squared,
         all_snp_r_squared = r_squared_snps, all_rand_r_squared = r_squared_rand) %>%
  filter(marker != "Q") # Remove population structure effect

## Copy the traits that were not included in the herit environment set
gwas_mlmm_model_adj <- bind_rows(
  gwas_mlmm_model_adj,
  filter(gwas_mlmm_model_adj, !trait %in% {distinct(gwas_mlmm_model_adj, set, trait) %>% split(.$set) %>% map("trait") %>% reduce(intersect)}) %>% 
    mutate(set = "herit_environments")
) %>%
  filter(set == "herit_environments")




# Iterate over each trait and term
pheno_fw_gen_corr_snp <- pheno_fw_use_tomodel %>%
  group_by(trait, term) %>%
  do({
    
    df <- .
    # Create a model.frame
    mf <- model.frame(estimate ~ line_name, df)
    
    # Create model matrices
    Z <- model.matrix(~ -1 + line_name, mf)
    Xmu <- model.matrix(~ 1, mf)
    
    ## Get the significant SNPs for this trait/term combination
    snps <- subset(gwas_mlmm_model_adj, trait == unique(df$trait) & coef %in% c("g", unique(df$term)), marker, drop = T)

    X_snp <- M[,snps,drop = F]
    X <- cbind(Xmu, X_snp)
    
    # Create the response matrix
    Y <- as.matrix(select(df, g, estimate))
    Y <- scale(Y)
    
    # Fit the model
    fit <- emmremlMultivariate(Y = t(Y), X = t(X), Z = t(Z), K = K)
    vcovG <- fit$Vg
    
    # Estimate the correlation using the formula for correlation
    # covariance / (sd1 * sd2)
    rhoG_hat <- vcovG[1,2] / prod(sqrt(diag(vcovG)))
    
    # Return
    data_frame(trait1 = "g", trait2 = unique(df$term), corG = rhoG_hat, df = nrow(df) - 2)
  })

# trait       term      trait1 trait2      corG    df
# 1 GrainYield  b         g      b          0.530   173
# 2 GrainYield  log_delta g      log_delta -0.735   173
# 3 HeadingDate b         g      b         -0.313   173
# 4 HeadingDate log_delta g      log_delta  0.225   173
# 5 PlantHeight b         g      b          0.380   173
# 6 PlantHeight log_delta g      log_delta  0.988   173



## Combine the phenotypic correlations and genetic correlations and output a table
stability_mean_corr_toprint <- pheno_fw_gen_corr_sig %>% 
  left_join(., stability_mean_corr_perm, by = c("trait", "term")) %>% 
  left_join(., pheno_fw_gen_corr_snp_sig, by = c("trait", "term", "trait1", "trait2")) %>%
  select(Trait = trait, Character1 = trait1, Character2 = trait2, corP = base, 
         pvalue_P = pvalue.y, corG = corG.x, pvalue_G = pvalue.x, corG_adj = corG.y, pvalue_Gadj = pvalue) %>%
  mutate_at(vars(contains("Character")), funs(str_replace_all(., coef_replace))) %>%
  mutate(PhenotypicCorrelation = str_c(round(corP, 3), " (p = ", formatC(pvalue_P, digits = 3, format = "e"), ")"),
         GeneticCorrelation = str_c(round(corG, 3), " (p = ", formatC(pvalue_G, digits = 3, format = "e"), ")"),
         AdjustedGeneticCorrelation = str_c(round(corG_adj, 3), " (p = ", formatC(pvalue_Gadj, digits = 3, format = "e"), ")")) %>%
  select(-contains("pvalue"), -contains("cor", ignore.case = FALSE))

# Print
write_csv(x = stability_mean_corr_toprint, path = file.path(fig_dir, "stability_mean_correlation.csv"))







# 
# # Count the number of cross-overs per genotype
# 
# # Group by trait
# pheno_fw_nXO <- pheno_fw_use_toplot %>% 
#   group_by(trait) %>% 
#   do({
#     df <- .
#     
#     # What is the range in environemntal values?
#     h_range <- range(df$h)
#     
#     # Get the distinct values of b and g for all lines
#     df1 <- distinct(df, line_name, g, b)
#     
#     # List of genotypes
#     line_names <- unique(df1$line_name)
#     
#     # Create a df
#     cross_over_count <- expand.grid(geno1 = line_names, geno2 = line_names) %>%
#       filter(geno1 != geno2)
#     
#     # Iterate over the data.frame
#     cross_over_count1 <- apply(X = cross_over_count, MARGIN = 1, FUN = function(i) {
#         
#       # Subset the information
#       df1_sub <- subset(df1, line_name %in% c(i), c(b, g))
#       
#       # Create a square matrix
#       A <- cbind(c(1, 1), -df1_sub$b)
#       b <- df1_sub$g
#       
#       # Solve
#       x <- solve(A, b)
#       
#       # Get the x value
#       x1 <- x[2]
#       
#       # Is the value within the environment range?
#       between(x1, h_range[1], h_range[2]) })
#     
#     # Combine the dfs
#     cbind(cross_over_count, cross_over_count1) })
# 
# 
# # Find the average proportion of crossovers per trait
# pheno_fw_nXO %>% 
#   group_by(trait, geno1) %>% 
#   summarize(nCO = sum(cross_over_count1), 
#             meanCO = mean(cross_over_count1)) %>% 
#   summarize(mean = mean(meanCO))
# 
# # trait        mean
# # 1 GrainYield  0.415
# # 2 HeadingDate 0.124
# # 3 PlantHeight 0.374
# 
# # Visualize
# pheno_fw_nXO %>% 
#   group_by(trait, geno1) %>% 
#   summarize(nCO = sum(cross_over_count1), 
#             meanCO = mean(cross_over_count1)) %>% 
#   ggplot(aes(x = trait, y = meanCO)) + 
#   geom_boxplot() + 
#   theme_bw()
# 






## Repeatability of stability using resampling

# Read in the results
load(file.path(result_dir, "pheno_fw_resampling.RData"))

# # Results from the TP
# pheno_sample_fw_tomodel <- pheno_samples_fw %>%
#   mutate(log_delta = log(delta)) %>%
#   select(-delta) %>%
#   gather(coef, value, b:log_delta)

# Results from TP + VP
pheno_sample_fw_tomodel <- pheno_samples_fw_tpvp %>%
  mutate(log_delta = log(delta)) %>%
  select(-delta) %>%
  gather(coef, value, g:log_delta)


## For each individual, calculate a coefficient of determination over all samples
estimate_cv <- pheno_sample_fw_tomodel %>% 
  mutate(p = as.factor(p)) %>%
  filter(!is.infinite(value),
         line_name != "07MT-10") %>% # Remove a line giving weird results
  group_by(trait, p, coef, line_name) %>% 
  summarize(estimate_cv = sd(value) / abs(mean(value)))

estimate_cv_summ <- estimate_cv %>%
  filter(estimate_cv <= 2) %>%
  summarize(mean_cv = mean(estimate_cv),
            lower = quantile(estimate_cv, alpha / 2),
            upper = quantile(estimate_cv, 1 - (alpha / 2))) # Calculate the mean variance and CI


## Plot mean and CI
g_estimate_cv_summ <- estimate_cv_summ %>%
  ungroup() %>%
  mutate(coef = str_replace_all(coef, coef_replace)) %>%
  ggplot(aes(x = p, y = mean_cv, ymin = lower, ymax = upper, group = 1)) + 
  geom_point(size = 0.2) +
  geom_line(stat = "identity", lwd = 0.2) +
  geom_errorbar(width = 0.3, lwd = 0.2) +
  ylab("Coefficient of variation") +
  xlab("Proportion of environments") + 
  facet_grid(coef~ trait, scales = "free_y", switch = "y") + 
  theme_pnas() +
  theme()

# Save
ggsave(filename = "pheno_fw_resample_cv.jpg", path = fig_dir, plot = g_estimate_cv_summ,
       height = 8, width = 8.7, units = "cm", dpi = 1000)


## Alternatively, find the mean squared deviation of each sample from the estimate using all environments
pheno_sample_fw_tomodel1 <- pheno_sample_fw_tomodel %>% 
  left_join(., distinct(pheno_mean_fw, trait, line_name, g, b, delta) %>%
              mutate(log_delta = log(delta)) %>% 
              select(-delta) %>% 
              gather(coef, true_value, g:log_delta))

# Calculate deviations, then summarize
estimate_deviations <- pheno_sample_fw_tomodel1 %>% 
  mutate(p = as.factor(p)) %>%
  filter(!is.infinite(value),
         line_name != "07MT-10") %>% # Remove a line giving weird results
  group_by(trait, p, coef, line_name) %>% 
  mutate(deviation = (value - true_value)^2) %>% 
  summarize(mean_deviation = mean(deviation))

estimate_deviations_summ <- estimate_deviations %>%
  summarize(mean_sqr_deviation = mean(mean_deviation),
            lower = quantile(mean_deviation, alpha / 2),
            upper = quantile(mean_deviation, 1 - (alpha / 2))) # Calculate the mean variance and CI

## Plot mean and CI
g_fw_deviation_summ <- estimate_deviations_summ %>%
  ungroup() %>%
  mutate(coef = str_replace_all(coef, coef_replace)) %>%
  ggplot(aes(x = p, y = mean_sqr_deviation, ymin = lower, ymax = upper, group = 1)) + 
  geom_point(size = 0.2) +
  geom_line(stat = "identity", lwd = 0.2) +
  geom_errorbar(width = 0.3, lwd = 0.2) +
  ylab("Mean squared deviation of estimate") +
  xlab("Proportion of environments") + 
  facet_wrap(coef~ trait, scales = "free_y", ncol = 3) + 
  theme_pnas() +
  theme()

# Save
ggsave(filename = "pheno_fw_resample_deviation.jpg", path = fig_dir, plot = g_fw_deviation_summ,
       height = 8, width = 8.7, units = "cm", dpi = 1000)












## Calculate the proportion of stability variance due to genomewide markers
# Load the results
load(file.path(result_dir, "stability_mean_marker_varcomp_results.RData"))

# Colors to use
color_use <- setNames(color[c(1, 3, 4)], c("Plastic", "Stable", "All Markers"))


## What is the proportion of variance explained by all genomewide markers
all_marker_prop_varcomp <- all_marker_varcomp %>% 
  ungroup() %>% 
  mutate(varP = varG + varR) %>% 
  mutate_at(vars(varG, varR), funs(. / varP))

all_marker_prop_varcomp %>%
  select(trait:varG) %>% 
  spread(coef, varG)

#   trait           b     g     log_delta
# <chr>       <dbl> <dbl>         <dbl>
# 1 GrainYield  0.436 0.437   0.0859 
# 2 HeadingDate 0.449 0.951   0.0275 
# 3 PlantHeight 0.289 0.749   0.00107

## Plastic and stable markers
plas_stab_marker_prop_varcomp <- plas_stab_marker_varcomp %>% 
  ungroup() %>% 
  mutate(varP = plastic + stable + residual) %>% 
  mutate_at(vars(-trait:-iter, -varP), funs(. / varP))

## Summarize
plas_stab_marker_prop_varcomp %>% 
  group_by(trait, coef) %>% 
  summarize_at(vars(plastic:residual), mean)

# trait       coef      plastic  stable residual
# 1 GrainYield  b          0.497  0          0.503
# 2 GrainYield  g          0.331  0.117      0.551
# 3 GrainYield  log_delta  0      0.0747     0.925
# 4 HeadingDate b          0.659  0          0.341
# 5 HeadingDate g          0.334  0.472      0.195
# 6 HeadingDate log_delta  0.0595 0.00104    0.939
# 7 PlantHeight b          0.462  0          0.538
# 8 PlantHeight g          0.0803 0.384      0.535
# 9 PlantHeight log_delta  0.0153 0.00364    0.981


## Plot
g_plas <- plas_stab_marker_prop_varcomp %>% 
  select(-residual, -varP) %>%
  gather(marker_type, varG, plastic, stable) %>%
  bind_rows(., mutate(all_marker_prop_varcomp, marker_type = "All Markers") %>% select(trait, coef, marker_type, varG)) %>%
  mutate(coef = str_replace_all(coef, coef_replace),
         marker_type = str_to_title(marker_type)) %>%
  ggplot(aes(x = coef, y = varG, color = marker_type)) + 
  geom_boxplot() +
  scale_color_manual(values = color_use, name = NULL) +
  facet_wrap(~ trait, ncol = 1, strip.position = "left") +
  ylab("Proportion of variance explained") +
  theme_pnas() +
  theme(legend.position = c(0.85, 0.54), axis.title.x = element_blank(), strip.placement = "outside",
        panel.spacing.y = unit(1, units = "line"))


# Save
ggsave(filename = "plastic_stable_maker_varcomp.jpg", plot = g_plas, path = fig_dir,
       height = 10, width = 8.7, units = "cm", dpi = 1000)










## Calculate marker effects and assess correlation
pheno_fw_mar_eff <- pheno_fw_use %>% 
  distinct(trait, line_name, g, term, estimate) %>% 
  spread(term, estimate) %>%
  gather(coef, value, g:log_delta) %>%
  group_by(trait, coef) %>% 
  do(data_frame(marker = colnames(M), effect = mixed.solve(y = .$value, Z = M)$u))

# Correlate
pheno_fw_mar_eff %>% 
  spread(coef, effect) %>% 
  summarize(cor_gb = cor(g, b),
            cor_gl = cor(g, log_delta),
            cor_bl = cor(b, log_delta))

# Plot
pheno_fw_mar_eff_toplot <- pheno_fw_mar_eff %>% 
  nest(marker, effect) %>% 
  crossing(., .) %>%
  filter(trait == trait1, coef != coef1) %>% 
  unnest()

g_marker_eff_corr <- pheno_fw_mar_eff_toplot %>%
  filter(coef == "g", coef1 != "g") %>%
  mutate_at(vars(contains("coef")), funs(str_replace_all(., coef_replace))) %>%
  ggplot(aes(x = effect, y = effect1)) +
  geom_point(size = 0.1) +
  facet_wrap(trait ~ coef1, scales = "free", ncol = 2) +
  theme_pnas() +
  theme(axis.title = element_blank())

## Save
ggsave(filename = "marker_effect_correlation.jpg", plot = g_marker_eff_corr, path = fig_dir,
       height = 8, width = 8.7, units = "cm", dpi = 1000)








#### Examples for presentatation


## Plot grain yield and HD
sample_traits <- c("GrainYield", "HeadingDate")

colors <- neyhart_palette("umn2")[c(3, 5)]
stab_cat <- c("Stable", "Average", "Plastic")


# Plot the normal FW analysis plot (genotype mean against environmental effects)
g_pheno_fw_b <- pheno_fw_use_toplot %>%
  filter(trait %in% sample_traits) %>%
  mutate(trait = str_replace_all(trait, traits_replace)) %>%
  split(.$trait) %>%
  map(~{
    df <- .
    
    # Create the breaks
    breaks <- pretty(range(df$b))
    breaks1 <- c(min(breaks), max(breaks))
    breaks2 <- c(1 - min(abs(breaks1 - 1)), 1, 1 + min(abs(breaks1 - 1)))
    limits <- c(1 - max(abs(breaks1 - 1)), 1, 1 + max(abs(breaks1 - 1)))
    
    y_breaks <- pretty(range(df$value))
    y_limits <- range(y_breaks)
    
    # Define colors
    cols <- colorRampPalette(c(colors[2], "grey75", colors[1]))(100)[c(1:45, 50, 55:100)]
    # Number of environments
    nE <- length(unique(df$environment))
    
    ggplot(data = df, aes(x = h, y = value, group = line_name)) + 
      geom_point(size = 1, color = "black") + 
      geom_abline(data = distinct(df, line_name, g, b), aes(slope = b, intercept = g, color = b), alpha = 0.8, lwd = 1) +
      facet_wrap(~ trait, nrow = 1, scales = "free") +
      scale_color_gradientn(colours = cols, name = NULL, breaks = breaks2, labels = stab_cat, limits = range(limits))  +
      ylab("Phenotypic value") +
      xlab("Environmental index") +
      scale_x_continuous(breaks = pretty) +
      scale_y_continuous(breaks = y_breaks, limits = y_limits) +
      theme_presentation() +
      theme(legend.position = c(0.20, 0.80), legend.direction = "vertical", legend.margin = margin(1, 1, 1, 1, "mm"),
            axis.title.x = element_blank())
    
  })

## Combine
g_fw_combine <- plot_grid(g_pheno_fw_b[[1]], g_pheno_fw_b[[2]] + theme(legend.position = "none", axis.title.y = element_blank()), 
                          nrow = 1, rel_widths = c(1, 0.85))
g_fw_combine1 <- ggdraw(add_sub(g_fw_combine, label = "Environmental index"))

ggsave(filename = "fw_plot_example.jpg", plot = g_fw_combine1, path = file.path(fig_dir, "Other"), width = 9, height = 5, dpi = 1000)


## Highlight examples of high and low stability
g_pheno_fw_b_example <- pheno_fw_use_toplot %>%
  filter(trait %in% sample_traits) %>%
  mutate(trait = str_replace_all(trait, traits_replace)) %>%
  split(.$trait) %>%
  map(~{
    df <- .
    
    # Create the breaks
    breaks <- pretty(range(df$b))
    breaks1 <- c(min(breaks), max(breaks))
    breaks2 <- c(1 - min(abs(breaks1 - 1)), 1, 1 + min(abs(breaks1 - 1)))
    limits <- c(1 - max(abs(breaks1 - 1)), 1, 1 + max(abs(breaks1 - 1)))
    
    y_breaks <- pretty(range(df$value))
    y_limits <- range(y_breaks)
    
    # Define colors
    cols <- colorRampPalette(c(colors[2], "grey75", colors[1]))(100)[c(1:45, 50, 55:100)]
    # Number of environments
    nE <- length(unique(df$environment))
    
    ## Subset the df for example genotypes
    n_sub <- 3
    df_distinct <- distinct(df, line_name, b)
    df1 <- bind_rows(top_n(df_distinct, n_sub, b), top_n(df_distinct, -n_sub, b)) %>%
      left_join(., df)
    
    
    ggplot(data = df1, aes(x = h, y = value, group = line_name)) + 
      geom_point(size = 1, color = "black") + 
      geom_abline(data = distinct(df1, line_name, g, b), aes(slope = b, intercept = g, color = b), alpha = 0.8, lwd = 1) +
      facet_wrap(~ trait, nrow = 1, scales = "free") +
      scale_color_gradientn(colours = cols, name = NULL, breaks = breaks2, labels = stab_cat, limits = range(limits))  +
      ylab("Phenotypic value") +
      xlab("Environmental index") +
      scale_x_continuous(breaks = pretty) +
      scale_y_continuous(breaks = y_breaks, limits = y_limits) +
      theme_presentation() +
      theme(legend.position = c(0.20, 0.80), legend.direction = "vertical", legend.margin = margin(1, 1, 1, 1, "mm"),
            axis.title.x = element_blank())
    
  })

## Combine
g_fw_example_combine <- plot_grid(g_pheno_fw_b_example[[1]], g_pheno_fw_b_example[[2]] + theme(legend.position = "none", axis.title.y = element_blank()), 
                                  nrow = 1, rel_widths = c(1, 0.85))
g_fw_example_combine1 <- ggdraw(add_sub(g_fw_example_combine, label = "Environmental index"))

ggsave(filename = "fw_plot_example_subset.jpg", plot = g_fw_example_combine1, path = file.path(fig_dir, "Other"), width = 9, height = 5, dpi = 1000)




colors <- setNames(object = neyhart_palette("umn1")[3:4], nm = coef_replace[1:2])

## Plot mean versus stability
g_mean_stability <- pheno_fw_use_toplot %>%
  distinct(trait, line_name, g, b, log_delta) %>%
  filter(trait %in% sample_traits) %>%
  mutate(trait = str_replace_all(trait, traits_replace)) %>%
  gather(term, estimate, b, log_delta) %>%
  mutate(term = str_replace_all(term, coef_replace)) %>%
  split(.$trait) %>%
  map(~{
    df <- .
    
    # Correlation
    df_cor <- df %>%
      group_by(term) %>% 
      do(test = cor.test(x = .$g, y = .$estimate)) %>% 
      ungroup() %>% 
      mutate(estimate = map_dbl(test, "estimate"), 
             p_value = map_dbl(test, "p.value"),
             p_value1 = ifelse(p_value < 0.05, "< 0.05", paste0("= ", formatC(p_value, 2))), 
             annotation = paste0("r = ", formatC(estimate, 2), " (P ", p_value1, ")"))
  
    ggplot(df, aes(x = g, y = estimate)) +
      geom_point(aes(color = term)) +
      geom_smooth(method = "lm", se = FALSE, color = "black") +
      # geom_label(data = df_cor, aes(x = Inf, y = -Inf, label = annotation), hjust = 1, vjust = -1, label.size = 0) + 
      facet_grid(term ~ trait, scales = "free_y", switch = "y") +
      ylab("Stability estimate") +
      xlab("Genotype mean") +
      scale_color_manual(values = colors, guide = FALSE) + 
      scale_y_continuous(breaks = pretty) +
      scale_x_continuous(breaks = pretty) +
      theme_presentation() +
      theme(strip.placement = "outside", axis.title.x = element_blank())
    
  })
    
    
## Combine
g_mean_stability_combine <- plot_grid(g_mean_stability[[1]], 
                                      g_mean_stability[[2]] + theme(axis.title.y = element_blank(), strip.text.y = element_blank()), 
                                      nrow = 1, rel_widths = c(1, 0.85))
g_mean_stability_combine1 <- ggdraw(add_sub(g_mean_stability_combine, label = "Genotype mean"))

ggsave(filename = "mean_stability_correlation_example.jpg", plot = g_mean_stability_combine1, path = file.path(fig_dir, "Other"), 
       width = 8, height = 5.5, dpi = 1000)









