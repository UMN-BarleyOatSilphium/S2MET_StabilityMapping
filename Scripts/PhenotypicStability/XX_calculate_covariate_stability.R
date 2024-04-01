## S2MET StabilityMapping
##
## Calculate the trait mean per se, overall stability, and stability for
## specific environmental covariates
##
##

# Load packages and set directories
library(modelr)
library(lme4)

# Project and other directories
source("startup.R")



# Load data ---------------------------------------------------------------

load(file.path(data_dir, "project_data.RData"))



# Calculate environmental means -------------------------------------------

# Convert factors
pheno_dat1 <- pheno_dat %>%
  mutate_at(vars(environment, line_name), as.factor)


## Fit the basic model with line name as random and environment as fixed.
## Because the environments are fixed, we don't really need to worry about
## specifying the R matrix correctly
## 
## Therefore, just use LMER


# Calculate environmental index --------------------------

# Iterate over traits
genotype_environment_mean <- pheno_dat1 %>%
  group_by(trait) %>%
  do({
    df <- .
    df1 <- droplevels(df) %>%
      mutate(environment = fct_contr_sum(environment))
    
    # Fit the model
    fit <- lmer(value ~ 1 + environment + (1|line_name), data = df1)
    
    calc_gh(object = fit)
    
  }) %>% ungroup()


# Calculate overall trait stability ---------------------------------------

# Regress the phenotypic value on the environmental mean
pheno_dat2 <- pheno_dat1 %>% 
  left_join(., genotype_environment_mean) %>%
  mutate(smith.weights = 1 / (std_error^2))

# Calculate linear and non-linear stability
genotype_stability <- pheno_dat2 %>%
  group_by(trait) %>%
  do(fwr(data = .)) %>%
  ungroup()



# Calculate stability in response to covariates ---------------------------

# Vector of covariate names
covariate_names <- names(ec_tomodel_centered$daymet)[-1:-2]

## Filter covariates depending on when the trait is observed
# Create a data.frame of covariates per trait
trait_covariate_df <- crossing(trait = traits, covariate = covariate_names) %>%
  filter(
    !(trait == "HeadingDate" & str_detect(covariate, "flowering|grain_fill")),
    !(trait == "PlantHeight" & str_detect(covariate, "grain_fill"))
  )


# First use a stepwise approach to add covariates
covariate_variable_selection_list <- genotype_environment_mean %>%
  left_join(., ec_tomodel_centered$daymet) %>%
  group_by(trait) %>%
  nest() %>%
  mutate(out = list(NULL))

# Iterate over rows
for (i in seq(nrow(covariate_variable_selection_list))) {

  tr <- covariate_variable_selection_list$trait[i]
  data <- covariate_variable_selection_list$data[[i]]

  # Subset trait-specific covariates
  covariates_use <- subset(trait_covariate_df, trait == tr, covariate, drop = TRUE)

  # Fit a base model
  fit_base <- lmer(value ~ (1|line_name) + (1|environment), data = data, weights = wts)

  ## Get the mean squares of the interaction
  int_meansq <- sigma(fit_base)^2

  # Add covariates one-at-a-time
  add1_list <- map(covariates_use,
                   ~update(fit_base, formula = add_predictors(formula(fit_base), reformulate(paste0("line_name:", .x)))))
  names(add1_list) <- covariates_use

  # Get the mean squares for those covariates
  add1_anova_list <- map(add1_list, ~anova(.x, type = "II"))

  # Find those covariates with significant interactions after multiple testing correction
  sig_ec_df <- add1_anova_list %>%
    map_df(tidy) %>%
    mutate(p_adj = p.adjust(p = p.value, method = "bonf")) %>%
    filter(p_adj < 0.01)

  sig_ecs <- str_remove(string = sig_ec_df$term, pattern = "line_name:")

  # Return a data.frame if there are significant ec
  if (length(sig_ecs) > 0) {

    sig_ec_betas <- add1_list[sig_ecs] %>%
      map(~fixef(.x) %>% tibble(term = names(.), estimate = .) %>%
            mutate(line_name = str_remove(map_chr(str_split(term, ":"), 1), "line_name")) %>%
            filter(str_detect(term, "Intercept", negate = TRUE)) %>% select(line_name, estimate)) %>%
      imap_dfr(~mutate(.x, covariate = paste0("response_", .y))) %>%
      spread(covariate, estimate) %>%
      left_join(distinct(data, line_name, g), .)

  } else {
    sig_ec_betas <- NA

  }

  ## Add to data frame and return
  covariate_variable_selection_list$out[[i]] <- sig_ec_betas

}


# Remove the traits with no covariates
covariate_variable_selection <- covariate_variable_selection_list %>%
  filter(sapply(out, is.data.frame)) %>%
  select(-data)





# Combine mean, stability, and covariate reaction data --------------------

stability_data_use <- inner_join(genotype_stability, covariate_variable_selection) %>%
  mutate(stability_info = map2(regression, out, inner_join, by = "line_name") %>%
           map(~select(., -genotype_mean) %>% rename(genotype_mean = g))) %>%
  select(trait, stability_info)




# Analyze stability -------------------------------------------------------


# Calculate genomic heritability of stability
stability_genomic_heritability <- stability_data_use %>%
  mutate(stability_info = map(stability_info, ~gather(.x, parameter, value, -line_name, -nEnv))) %>%
  unnest() %>%
  group_by(trait, parameter) %>%
  do({
    df <- .

    # Calculate a residual matrix
    # If we are looking at the genotype means, use the reciprocal of the number
    # of environments
    # For stability, use 1.
    # if (unique(df$parameter) == "genotype_mean") r <- df$nEnv else r <- 1

    R <- diag(1, nrow = nrow(K))
    dimnames(R) <- dimnames(K)

    # Subset data
    df1 <- filter(df, line_name %in% row.names(K))

    # Calculate genomic heritability
    std_out <- capture.output({
      genomic_h2 <- marker_h2_means(data.vector = df1$value, geno.vector = df1$line_name,
                                    K = K, Dm = R, alpha = 0.05, max.iter = 1000)
    })

    ## Return heritability with CI
    tibble(h2 = genomic_h2$h2, lower = genomic_h2$conf.int1[1], upper = genomic_h2$conf.int1[2])

  }) %>% ungroup()



# Save
save("stability_data_use", "stability_genomic_heritability", "genotype_environment_mean", "genotype_stability",
     file = file.path(result_dir, "stability_estimates.RData"))

