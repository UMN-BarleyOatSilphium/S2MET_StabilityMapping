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
library(parallel)

# Repository directory
repo_dir <- getwd()
# Project and other directories
source(file.path(repo_dir, "source.R"))
 
# Load fw data
load(file.path(result_dir, "pheno_mean_fw_results.RData"))

# Relationship matrix
M <- S2TP_imputed_multi_genos_mat
K <- A.mat(X = M, min.MAF = 0, max.missing = 1)

# Significance threshold
alpha <- 0.05

# Number of cores
n_cores <- detectCores()

# Substitute the estimates of stability using just the TP with estimates using
# both the TP and VP
pheno_mean_fw <- pheno_mean_fw_tpvp %>%
  filter(line_name %in% tp)


### Preliminary analysis

## Are traits correlated with latitude?
env_means <- pheno_mean_fw %>% 
  distinct(trait, environment, h) %>%
  left_join(distinct(trial_info, environment, latitude, longitude))

env_means %>%
  group_by(trait) %>%
  summarize_at(vars(latitude, longitude), funs(corr = list(neyhart::bootstrap(x = h, y = ., fun = "cor")))) %>%
  gather(coord, out, -trait) %>%
  unnest() %>% 
  mutate(significant = ! (ci_lower <= 0 & ci_upper >= 0))





### Finlay-Wilkinson Regression

## Results

# What are the ranges of genotype values and environmental effects per trait?
pheno_summary <- pheno_mean_fw %>% 
  distinct(trait, line_name, environment, h, g)  %>% 
  gather(term, value, g:h) %>% 
  group_by(trait, term) %>% 
  summarize_at(vars(value), funs(min, max)) %>% 
  mutate(range = max - min)

## Environment metadata
trial_info %>% 
  filter(environment %in% pheno_mean_fw$environment) %>%
  select(environment, latitude, longitude) %>%
  gather(measure, value, -environment) %>% 
  group_by(measure) %>% 
  summarize_at(vars(value), funs(min, max), na.rm = T) %>% 
  mutate(range = max - min)


### Phenotypic Stability

# Transformation

# Extract the unique stability coefficients
# Then add the breeding program information
pheno_mean_fw1 <- pheno_mean_fw %>%
  left_join(., subset(entry_list, Class == "S2TP", c(Line, Program)), 
            by = c("line_name" = "Line")) %>%
  rename(program = Program)

## Log transform the non-linear stability estimates
pheno_fw_use <- pheno_mean_fw1 %>%
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
  theme_bw() +
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
g_fw_dist <- plot_grid(g_pheno_fw_dens_b, g_pheno_fw_dens_delta, ncol = 2)

ggsave(filename = "stability_estimate_distriubtions.jpg", plot = g_fw_dist, path = fig_dir,
       height = 10, width = 5, dpi = 1000)




# Plot the stability estimates as lines against g_i vs h_j
# colors_use <- set_names(umn_palette(2)[3:7], unique(pheno_fw_use$program))
color_use <- umn_palette(2)[3]

# Add the program information to the results
pheno_fw_use_toplot <- pheno_fw_use %>% 
  spread(term, estimate) %>% 
  select(trait, environment, line_name, value, g, h, b) %>% 
  left_join(., subset(entry_list, Class == "S2TP", c(Line, Program)), by = c("line_name" = "Line")) %>%
  rename(program = Program)

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



## Plot the relationship between the genotypic effect and the sensitivity
set.seed(415)

# Calculate the correlation between genotype mean and the stability estimate
stability_mean_corr <- pheno_fw_use %>%
  distinct(trait, line_name, g, term, estimate) %>%
  group_by(trait, term) %>% 
  # Bootstrap the correlation
  do(neyhart::bootstrap(x = .$g, y = .$estimate, fun = "cor")) %>%
  mutate(significant = !between(0, ci_lower, ci_upper),
         sig_ann = ifelse(significant, "*", ""),
         annotation = str_c("r = ", round(base, 3), sig_ann))


# Use a parametric test
stability_mean_corr_param <- pheno_fw_use %>%
  distinct(trait, line_name, g, term, estimate) %>%
  group_by(trait, term) %>% do(test = cor.test(.$g, .$estimate)) %>% 
  ungroup() %>% 
  mutate(base = map_dbl(test, "estimate"), 
         pvalue = map_dbl(test, "p.value")) %>% 
  select(-test) %>%
  mutate(significant = pvalue <= alpha,
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
  filter(term == "b") %>%
  ggplot(aes(x = g, y = estimate, group = FALSE)) +
  geom_text(data = subset(stability_mean_corr_param, term == "b"), 
            aes(x = Inf, y = -Inf, label = annotation, vjust = -1, hjust = 1.2), col = "black", size = 1.5) +
  facet_wrap(~ trait, scales = "free_x", ncol = 1, strip.position = "left") +
  ylab("Linear stability estimate") +
  theme(legend.position = "right") +
  g_add

# Plot just the non-linear stability
g_nonlinear_stability_and_mean <- pheno_fw_use %>%
  filter(term == "log_delta") %>%
  ggplot(aes(x = g, y = estimate, group = FALSE)) +
  geom_text(data = subset(stability_mean_corr_param, term == "log_delta"),
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
library(EMMREML)

# Create a df for modeling
pheno_fw_use_tomodel <- pheno_fw_use %>%
  ungroup() %>%
  distinct(trait, line_name, g, term, estimate)

# Iterate over each trait and term
pheno_fw_gen_corr <- pheno_fw_use_tomodel %>%
  group_by(trait, term) %>%
  do({
    df <- .
    # Create a model.frame
    mf <- model.frame(estimate ~ line_name, df)
    
    # Create model matrices
    Z <- model.matrix(~ -1 + line_name, mf)
    X <- model.matrix(~ 1, mf)
    
    # Create the response matrix
    Y <- as.matrix(select(df, g, estimate))
    
    # Fit the model
    fit <- emmremlMultivariate(Y = t(Y), X = t(X), Z = t(Z), K = K)
    vcovG <- fit$Vg
    
    # Estimate the correlation using the formula for correlation
    # covariance / (sd1 * sd2)
    rhoG_hat <- vcovG[1,2] / prod(sqrt(diag(vcovG)))
    
    # Return
    data_frame(trait1 = "g", trait2 = unique(df$term), corG = rhoG_hat, df = nrow(df) - 2)
  })


## Use a t-distribution to compute a p-value
pheno_fw_gen_corr_sig <- pheno_fw_gen_corr %>% 
  ungroup() %>% 
  mutate(stat = sqrt(df) * corG/sqrt(1 - corG^2), 
         pvalue = 2 * pt(q = abs(stat), df = df, lower.tail = FALSE))

# trait       term      trait1 trait2       corG    df    stat   pvalue
# 1 GrainYield  b         g      b          0.646    181  11.4   5.13e-23
# 2 GrainYield  log_delta g      log_delta  0.310    181   4.38  1.99e- 5
# 3 HeadingDate b         g      b          0.0217   181   0.291 7.71e- 1
# 4 HeadingDate log_delta g      log_delta -0.289    181  -4.07  7.10e- 5
# 5 PlantHeight b         g      b          0.472    181   7.19  1.61e-11
# 6 PlantHeight log_delta g      log_delta  0.566    181   9.23  7.30e-17


## Combine the phenotypic correlations and genetic correlations and output a table
stability_mean_corr_toprint <- pheno_fw_gen_corr_sig %>% 
  left_join(., stability_mean_corr_param, by = c("trait", "term")) %>% 
  select(Trait = trait, Character1 = trait1, Character2 = trait2, corP = base, 
         pvalue_P = pvalue.y, corG, pvalue_G = pvalue.x) %>%
  mutate_at(vars(contains("Character")), funs(str_replace_all(., coef_replace))) %>%
  mutate(PhenotypicCorrelation = str_c(round(corP, 3), " (p = ", formatC(pvalue_P, digits = 3, format = "e"), ")"),
         GeneticCorrelation = str_c(round(corG, 3), " (p = ", formatC(pvalue_G, digits = 3, format = "e"), ")")) %>%
  select(-contains("pvalue"), -contains("cor", ignore.case = FALSE))

# Print
write_csv(x = stability_mean_corr_toprint, path = file.path(fig_dir, "stability_mean_correlation.csv"))


## Calculate the genetic correlation after correcting for significant GWAS markers
# Load the GWAS results
load(file.path(result_dir, "pheno_fw_mean_gwas_results.RData"))



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
    # snps <- subset(gwas_mlmm_final_tpvp, model == "QG" & trait == unique(df$trait) & coef == unique(df$term), term, drop = T)
    snps <- subset(gwas_mlmm_final_tpvp, model == "QG" & trait == unique(df$trait) & coef %in% c("g", unique(df$term)), term, drop = T)
    
    X_snp <- M[,snps]
    X <- cbind(Xmu, X_snp)
    
    # Create the response matrix
    Y <- as.matrix(select(df, g, estimate))
    
    # Fit the model
    fit <- emmremlMultivariate(Y = t(Y), X = t(X), Z = t(Z), K = K)
    vcovG <- fit$Vg
    
    # Estimate the correlation using the formula for correlation
    # covariance / (sd1 * sd2)
    rhoG_hat <- vcovG[1,2] / prod(sqrt(diag(vcovG)))
    
    # Return
    data_frame(trait1 = "g", trait2 = unique(df$term), corG = rhoG_hat, df = nrow(df) - 2)
  })


## Use a t-distribution to compute a p-value
pheno_fw_gen_corr_snp_sig <- pheno_fw_gen_corr_snp %>% 
  ungroup() %>% 
  mutate(stat = sqrt(df) * corG/sqrt(1 - corG^2), 
         pvalue = 2 * pt(q = abs(stat), df = df, lower.tail = FALSE))

# trait       term      trait1 trait2      corG    df   stat   pvalue
# 1 GrainYield  b         g      b          0.646   181  11.4  5.13e-23
# 2 GrainYield  log_delta g      log_delta  0.310   181   4.38 1.99e- 5
# 3 HeadingDate b         g      b         -0.443   181  -6.66 3.24e-10
# 4 HeadingDate log_delta g      log_delta -0.426   181  -6.34 1.80e- 9
# 5 PlantHeight b         g      b          0.412   181   6.08 6.96e- 9
# 6 PlantHeight log_delta g      log_delta  0.700   181  13.2  2.87e-28





# Count the number of cross-overs per genotype

# Group by trait
pheno_fw_nXO <- pheno_fw_use_toplot %>% 
  group_by(trait) %>% 
  do({
    df <- .
    
    # What is the range in environemntal values?
    h_range <- range(df$h)
    
    # Get the distinct values of b and g for all lines
    df1 <- distinct(df, line_name, g, b)
    
    # List of genotypes
    line_names <- unique(df1$line_name)
    
    # Create a df
    cross_over_count <- expand.grid(geno1 = line_names, geno2 = line_names) %>%
      filter(geno1 != geno2)
    
    # Iterate over the data.frame
    cross_over_count1 <- apply(X = cross_over_count, MARGIN = 1, FUN = function(i) {
        
      # Subset the information
      df1_sub <- subset(df1, line_name %in% c(i), c(b, g))
      
      # Create a square matrix
      A <- cbind(c(1, 1), -df1_sub$b)
      b <- df1_sub$g
      
      # Solve
      x <- solve(A, b)
      
      # Get the x value
      x1 <- x[2]
      
      # Is the value within the environment range?
      between(x1, h_range[1], h_range[2]) })
    
    # Combine the dfs
    cbind(cross_over_count, cross_over_count1) })


# Find the average proportion of crossovers per trait
pheno_fw_nXO %>% 
  group_by(trait, geno1) %>% 
  summarize(nCO = sum(cross_over_count1), 
            meanCO = mean(cross_over_count1)) %>% 
  summarize(mean = mean(meanCO))

# trait        mean
# 1 GrainYield  0.415
# 2 HeadingDate 0.124
# 3 PlantHeight 0.374

# Visualize
pheno_fw_nXO %>% 
  group_by(trait, geno1) %>% 
  summarize(nCO = sum(cross_over_count1), 
            meanCO = mean(cross_over_count1)) %>% 
  ggplot(aes(x = trait, y = meanCO)) + 
  geom_boxplot() + 
  theme_bw()







## Repeatability of stability using resampling

# Read in the results
load(file.path(result_dir, "pheno_fw_resampling.RData"))

pheno_sample_fw_tomodel <- pheno_sample_mean_fw %>%
  mutate(log_delta = log(delta)) %>%
  select(-delta) %>%
  gather(coef, value, b:log_delta)

# Fit a fixed-effect model
repeatability_fit <- pheno_sample_fw_tomodel %>% 
  filter(!is.infinite(value)) %>% 
  # filter(iter %in% c(1:40)) %>%
  group_by(trait, p, coef) %>% 
  do(fit = lm(value ~ line_name, .))

## Summarize
repeatability_summ <- repeatability_fit %>% 
  ungroup() %>% 
  mutate(anova = map(fit, ~tidy(anova(.)))) %>% 
  unnest(anova) %>% 
  group_by(trait, p, coef) %>% 
  summarize(repeatability = meansq[1] / sum(meansq))

# Plot
g_fw_repeatability <- repeatability_summ %>% 
  ggplot(aes(x = p, y = repeatability, group = 1)) + 
  geom_point() + 
  geom_line(stat = "identity") + 
  ylab("Repeatability") +
  xlab("Proportion of Environments") + 
  facet_grid(coef~ trait) + 
  theme_bw()

# Save
ggsave(filename = "pheno_fw_resample_repeatability.jpg", path = fig_dir, plot = g_fw_repeatability,
       height = 4, width = 6, dpi = 1000)



## Calculate the proportion of stability variance due to genomewide markers
# Load the results
load(file.path(result_dir, "stability_mean_marker_varcomp_results.RData"))

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
# 1 GrainYield  0.431 0.345 0.0813       
# 2 HeadingDate 0.449 0.611 0.0774       
# 3 PlantHeight 0.289 0.517 0.00000000100

## Random markers
rand_marker_prop_varcomp <- rand_marker_varcomp %>%
  mutate(varP = varG + varR,
         nmar = parse_number(nmar)) %>% 
  mutate_at(vars(varG, varR), funs(. / varP)) %>%
  group_by(trait, coef, nmar) %>%
  summarize_at(vars(varG), funs(mean = mean(.), lower = quantile(., alpha / 2), upper = quantile(., 1 - (alpha / 2)))) %>%
  ungroup()

## Plot
g_rand <- rand_marker_prop_varcomp %>% 
  ggplot(aes(x = nmar, y = mean, color = coef)) + 
  geom_point(pch = 0) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.1) + 
  facet_wrap(~trait, ncol = 2) +
  theme_bw()


## Calculate marker effects and assess accuracy
pheno_fw_me <- pheno_fw_use %>% 
  distinct(trait, line_name, term, estimate) %>% 
  filter(term == "b") %>%
  group_by(trait) %>% 
  do(me = mixed.solve(y = .$estimate, Z = M)$u)

# Predict






## Compare different marker subsets for the proportion of variance explained
# Top markers
top_marker_prop_varcomp <- top_marker_varcomp %>%
  mutate(varP = varG + varR,
         nmar = parse_number(nmar)) %>% 
  mutate_at(vars(varG, varR), funs(. / varP))

## Plot
g_top <- top_marker_prop_varcomp %>% 
  ggplot(aes(x = nmar, y = varG, color = coef)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~trait, ncol = 3) +
  theme_bw()

# Add random
g_top1 <- g_top + 
  geom_point(data = rand_marker_prop_varcomp, aes(y = mean), pch = 0) +
  geom_line(data = rand_marker_prop_varcomp, aes(y = mean)) +
  geom_ribbon(data = rand_marker_prop_varcomp, aes(ymin = lower, ymax = upper, y = mean), alpha = 0.1) +
  theme(legend.position = c(0.25, 0.80)) + 
  ylim(c(0, 1))
  
# Save
ggsave(filename = "top_maker_varcomp.jpg", plot = g_top1, path = fig_dir,
       height = 4, width = 8, dpi = 1000)
  
  

# evenly-spaced markers
esm_marker_prop_varcomp <- esm_marker_varcomp %>%
  mutate(varP = varG + varR,
         nmar = parse_number(nmar)) %>% 
  mutate_at(vars(varG, varR), funs(. / varP))

## Plot
g_esm <- esm_marker_prop_varcomp %>% 
  ggplot(aes(x = nmar, y = varG, color = coef)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~trait, ncol = 3) +
  theme_bw()

g_esm1 <- g_esm + 
  geom_point(data = rand_marker_prop_varcomp, aes(y = mean), pch = 0) +
  geom_line(data = rand_marker_prop_varcomp, aes(y = mean)) +
  geom_ribbon(data = rand_marker_prop_varcomp, aes(ymin = lower, ymax = upper, y = mean), alpha = 0.1) +
  theme(legend.position = c(0.25, 0.80)) + 
  ylim(c(0, 1))

# Save
ggsave(filename = "even_spaced_maker_varcomp.jpg", plot = g_esm1, path = fig_dir,
       height = 4, width = 8, dpi = 1000)





# Top evenly-spaced markers
tesm_marker_prop_varcomp <- tesm_marker_varcomp %>%
  mutate(varP = varG + varR,
         nmar = parse_number(nmar)) %>% 
  mutate_at(vars(varG, varR), funs(. / varP))

## Plot
g_tesm <- tesm_marker_prop_varcomp %>% 
  ggplot(aes(x = nmar, y = varG, color = coef)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~trait, ncol = 3) +
  theme_bw()

g_tesm1 <- g_tesm + 
  geom_point(data = rand_marker_prop_varcomp, aes(y = mean), pch = 0) +
  geom_line(data = rand_marker_prop_varcomp, aes(y = mean)) +
  geom_ribbon(data = rand_marker_prop_varcomp, aes(ymin = lower, ymax = upper, y = mean), alpha = 0.1) +
  theme(legend.position = c(0.25, 0.80)) + 
  ylim(c(0, 1))

# Save
ggsave(filename = "top_even_spaced_maker_varcomp.jpg", plot = g_tesm1, path = fig_dir,
       height = 4, width = 8, dpi = 1000)


# Top plastic markers
plas_marker_prop_varcomp <- plas_marker_varcomp %>%
  mutate(varP = varG + varR,
         nmar = parse_number(nmar)) %>% 
  mutate_at(vars(varG, varR), funs(. / varP))

## Plot
g_plas <- plas_marker_prop_varcomp %>% 
  ggplot(aes(x = nmar, y = varG, color = coef)) + 
  geom_point() + 
  geom_line() + 
  facet_wrap(~trait, ncol = 3) +
  theme_bw()

g_plas1 <- g_plas + 
  geom_point(data = rand_marker_prop_varcomp, aes(y = mean), pch = 0) +
  geom_line(data = rand_marker_prop_varcomp, aes(y = mean)) +
  geom_ribbon(data = rand_marker_prop_varcomp, aes(ymin = lower, ymax = upper, y = mean), alpha = 0.1) +
  theme(legend.position = c(0.25, 0.80)) + 
  ylim(c(0, 1))

# Save
ggsave(filename = "plastic_maker_varcomp.jpg", plot = g_plas1, path = fig_dir,
       height = 4, width = 8, dpi = 1000)






