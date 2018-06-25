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

# Repository directory
repo_dir <- getwd()
# Project and other directories
source(file.path(repo_dir, "source.R"))
 
# Load fw data
load(file.path(result_dir, "pheno_mean_fw_results.RData"))


# Significance threshold
alpha <- 0.05

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
  dplyr::rename(program = Program)

## Log transform the non-linear stability estimates
pheno_mean_fw_trans <- pheno_mean_fw1 %>%
  group_by(trait) %>% 
  mutate(log_delta = log(delta)) %>%
  # Tidy
  gather(term, estimate, b, delta, log_delta) %>% 
  filter(term != "delta")

# # Remove potential outliers (visually)
# pheno_mean_fw_trans_filter <- pheno_mean_fw_trans %>%
#   # Remove 07MT-10 from grain protein linear stability
#   filter(!(trait == "GrainProtein" & line_name == "07MT-10" & term == "b")) %>%
#   # Remove 07AB-84 and 08AB-08 from heading date non-linear stability
#   filter(!(trait == "HeadingDate" & term == "log_delta" & estimate > 3)) %>%
#   # Remove 07MT-10 from TestWeight
#   filter(!(trait == "TestWeight" & line_name == "07MT-10"))
  
pheno_fw_use <- pheno_mean_fw_trans %>%
  filter(trait %in% c("GrainYield", "HeadingDate", "PlantHeight"))

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

# Add the program information to the results
pheno_fw_use_toplot <- pheno_fw_use %>% 
  spread(term, estimate) %>% 
  select(trait, environment, line_name, value, g, h, b) %>% 
  left_join(., subset(entry_list, Class == "S2TP", c(Line, Program)), by = c("line_name" = "Line")) %>%
  dplyr::rename(program = Program)

# Plot the normal FW analysis plot (genotype mean against environmental effects)
g_pheno_fw_b <- pheno_fw_use_toplot %>%
  ggplot(aes(x = h, y = value, col = program, group = line_name)) + 
  geom_point() + 
  geom_abline(aes(slope = b, intercept = g, col = program), alpha = 0.15) +
  facet_wrap(~ trait, ncol = 1, scales = "free") +
  scale_color_discrete(drop = FALSE, guide = guide_legend(title = "Program", ncol = 1)) +
  ylab("Phenotypic Value") +
  xlab("Environmental Effect") +
  labs(title = "Linear Phenotypic Stability") +
  # scale_color_brewer(palette = "Set2") +
  theme_bw() +
  theme(legend.position = "right",
        legend.background = element_rect(fill = "grey85", color = NA))

ggsave(filename = "pheno_fw_linear.jpg", plot = g_pheno_fw_b, path = fig_dir, 
       width = 5, height = 10, dpi = 1000)


## Subset genotypes to highlight different responses
pheno_fw_example <- pheno_fw_use_toplot %>%
  distinct(trait, line_name, g, b) %>%
  group_by(trait) %>% 
  filter(b == min(b) | b == max(b) | abs(1 - b) == min(abs(1 - b)))

g_pheno_fw_b_example <- pheno_fw_use_toplot %>%
  select(trait, environment, line_name, value, h) %>% 
  ggplot(aes(x = h, y = value)) + 
  geom_point() + 
  geom_abline(data = pheno_fw_example, aes(intercept = g, slope = b, col = line_name), lwd = 1) + 
  scale_color_discrete(guide = guide_legend(title = "Line Name", ncol = 1)) +
  facet_wrap(~ trait, ncol = 1, scales = "free") +
  ylab("Phenotypic Value") +
  xlab("Environment Effect") +
  labs(title = "Example Genotype Responses") +
  theme_bw() +
  theme(legend.position = "right",
        legend.background = element_rect(fill = "grey85", color = NA),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10))

ggsave(filename = "pheno_fw_linear_example.jpg", plot = g_pheno_fw_b_example, path = fig_dir, 
       width = 5, height = 10, dpi = 1000)


# Combine plots
g_pheno_fw_figure <- plot_grid(g_pheno_fw_b, g_pheno_fw_b_example, labels = c("A", "B"), ncol = 2)

ggsave(filename = "pheno_fw_linear_combined.jpg", plot = g_pheno_fw_figure, path = fig_dir,
       width = 10, height = 10, dpi = 1000)



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

# Visualize
pheno_fw_nXO %>% 
  group_by(trait, geno1) %>% 
  summarize(nCO = sum(cross_over_count1), 
            meanCO = mean(cross_over_count1)) %>% 
  ggplot(aes(x = trait, y = meanCO)) + 
  geom_boxplot() + 
  theme_bw()




## Plot the relationship between the genotypic effect and the sensitivity
set.seed(415)

# Calculate the correlation between genotype mean and the stability estimate
stability_mean_corr <- pheno_fw_use %>%
  group_by(trait, term) %>% 
  # Bootstrap the correlation
  do(boot = {
    df <- .
    # Perform bootstrap
    boot_out <- df %>% 
      bootstrap(n = 1000) %>% 
      pull(strap) %>%
      map_dbl(~as.data.frame(.) %>% ungroup() %>% select(g, estimate) %>% cor(.) %>% .[1,2])
    
    # Summarize and export
    data.frame(cor = boot_out) %>% 
      summarize(boot_mean = mean(cor), se = sd(cor), 
                ci_lower = quantile(cor, alpha / 2), ci_upper = quantile(cor, 1 - (alpha / 2))) }) %>%
  # Join with original data
  left_join(pheno_fw_use, .) %>%
  # Calculate the base correlation and add a character annotation
  group_by(trait, term) %>% 
  mutate(corr = cor(g, estimate, use = "complete.obs"),
         corr = str_c("r = ", round(corr, 3))) %>% 
  unnest() %>%
  ungroup()


# # Fit linear models
# stability_mean_lm <- stability_mean_corr %>% 
#   group_by(trait, term) %>%
#   do(fit = lm(g ~ estimate, data = .))
# 
# # Tidy up
# stability_mean_lm_tidy <- stability_mean_lm %>%
#   ungroup() %>%
#   mutate(tidy_fit = map(fit, tidy)) %>%
#   unnest(tidy_fit) %>% 
#   filter(term1 == "estimate")


## Looks like significant regression slopes for the relationship between genotype 
## mean and $b$ for GrainYield and PlantHeight, and between genotype mean and $\delta$ 
## for HeadingDate and PlantHeight.

# The relationship between genotype mean and $\delta$ is still significant after outlier removal

# Extract the data to plot
stability_mean_corr_toplot <- stability_mean_corr %>% 
  rowwise() %>%
  mutate(significant = !between(0, ci_lower, ci_upper)) %>% 
  ungroup() %>%
  distinct(trait, line_name, program, g, term, estimate, corr, significant) %>%
  mutate(annotation = str_c(corr, ifelse(significant, "*", "")))

# Create a list of plot additions
g_add <- list(geom_point(),
              geom_smooth(method = "lm", se = FALSE, col = "black"),
              geom_text(aes(x = Inf, y = -Inf, label = annotation, vjust = -1, hjust = 1.2), col = "black"),
              scale_color_discrete(name = "Program"),
              ylab("Estimate"),
              xlab("Genotype Mean"),
              theme_bw(),
              theme(legend.background = element_rect(fill = "grey85", color = NA)))

# Plot just the linear stability
g_linear_stability_and_mean <- stability_mean_corr_toplot %>%
  filter(term == "b") %>%
  ggplot(aes(x = g, y = estimate, col = program, group = FALSE)) +
  facet_wrap(~ trait, scales = "free_x", ncol = 1) +
  labs(title = "Linear Stability") +
  theme(legend.position = "right") +
  g_add
  
# Plot just the non-linear stability
g_nonlinear_stability_and_mean <- stability_mean_corr_toplot %>%
  filter(term == "log_delta") %>%
  ggplot(aes(x = g, y = estimate, col = program, group = FALSE)) +
  facet_wrap(~ trait, scales = "free", ncol = 1) +
  labs(title = "Non-Linear Stability" ) +
  theme(legend.position = "right") +
  g_add
  

# Plot grid
# Remove the legend from each plot object. Also remove the y axis title from the right-hand plot
g_plotgrid <- plot_grid(
  g_linear_stability_and_mean + theme(legend.position="none", axis.title.x = element_blank()),
  g_nonlinear_stability_and_mean + theme(legend.position="none", axis.title.x = element_blank()) + ylab(""),
  align = "h",
  labels = c("A", "B"))

# Extract the legend from one of the plots and add it back in
g_plotgrid1 <- plot_grid(
  g_plotgrid, get_legend(g_linear_stability_and_mean), rel_widths = c(2, 0.3)
)

# Add a common x-axis legend.
g_plotgrid2 <- ggdraw(add_sub(g_plotgrid1, "Genotype Mean", size = 12, hjust = 0.7))

ggsave(filename = "pheno_fw_stability_mean.jpg", plot = g_plotgrid2, path = fig_dir, 
       width = 8, height = 6, dpi = 1000)



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
load(file.path(result_dir, "pheno_mar_fw_varcomp.RData"))

prop_varcomp <- pheno_fw_var_comp_all_markers %>% 
  mutate(total = snps + res) %>% 
  mutate_at(vars(snps, res), funs(. / total))

prop_varcomp %>%
  select(trait:snps) %>% 
  spread(coef, snps)

#   trait           b     g     log_delta
# <chr>       <dbl> <dbl>         <dbl>
# 1 GrainYield  0.431 0.345 0.0813       
# 2 HeadingDate 0.449 0.611 0.0774       
# 3 PlantHeight 0.289 0.517 0.00000000100



