## S2MET Mapping
## Calculate marker effects in each environment
## 
## This script will calculate genomewide marker effects in each environment, then
## calculate marker effect stability
## 
## Author: Jeff Neyhart
## Last updated: July 11, 2018
## 

# Run the source script - local
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))


# Rename the marker matrix
M <- s2tp_genos_imputed

# Get the marker information
snp_info <- s2tp_genos_hmp %>%
  select(marker, chrom:cM_pos) %>%
  # Correct the cM position for BOPA snps
  mutate(cM_pos = if_else(str_detect(marker, "^S"), cM_pos, cM_pos / 1000))


# Subset the entry list for the tp
tp_entry_list <- entry_list %>%
  filter(Line %in% tp_geno) %>%
  select(line_name = Line, program = Program)

# Significance threshold for GWAS
alpha <- 0.05

# What model will we use for mapping?
model_use <- "QG"


## Load the phenotypic stability results
load(file.path(result_dir, "pheno_mean_fw_results.RData"))

# Use the environment mean estimates from the TP and VP
pheno_mean_fw <- pheno_mean_fw_tpvp %>% 
  filter(line_name %in% tp_geno)


## Calculate marker effects separately for each environment
marker_by_env_effects <- pheno_mean_fw %>%
  mutate(line_name = as.factor(line_name)) %>%
  group_by(trait, environment) %>%
  do({
    df <- .

    # Model frame
    mf <- model.frame(value ~ line_name + std_error, df)

    # Response
    y <- model.response(mf)
    # Grand mean design
    X <- model.matrix(~ 1, mf)
    # Line_name design
    Zg <- model.matrix(~ -1 + line_name, mf)

    # Subset the marker matrix
    Z <- Zg %*% M
    K <- diag(ncol(Z))

    # Pull out the weights into an R matrix
    # R <- solve(diag(mf$std_error^2))

    # Fit
    # fit <- sommer::mmer(Y = y, X = X, Z = list(marker = list(Z = Z, K = K)), R = list(res = R))
    fit <- mixed.solve(y = y, Z = Z, method = "REML")


    # # Grab the marker effects and return
    # marker_effect <- fit$u.hat$marker %>%
    #   as.data.frame() %>%
    #   rownames_to_column("marker") %>%
    #   dplyr::select(marker, effect = T1)

    # Grab the marker effects and return
    marker_effect <- fit$u %>%
      as.data.frame() %>%
      rownames_to_column("marker") %>%
      rename(effect = ".")

    # Return a data_frame
    data_frame(marker_effect = list(marker_effect), fit = list(fit)) })




## Manage the mxe data and convert to a tidy data.frame
mxe_df <- marker_by_env_effects %>%
  unnest(marker_effect) %>%
  ungroup()

# Combine with the environmental mean
mxe_df1 <- mxe_df %>%
  left_join(., distinct(pheno_mean_fw, trait, environment, h),
            by = c("environment", "trait"))

# Calculate the marker effect stability
# Regress the marker effects in each environment on the mean in that environment
marker_effect_stab <- mxe_df1 %>%
  group_by(trait, marker) %>%
  do({
    df <- .
    # Fit the model
    fit <- lm(effect ~ h, df)

    # Output a data.frame with the intercept, slope, MSE, and fitted marker effects
    mutate(df, a = coef(fit)[1], c = coef(fit)[2], delta = mean(resid(fit)^2),
           fitted_effects = fitted(fit))

  })


# Add the snp information
marker_mean_fw <- marker_effect_stab %>%
  left_join(., snp_info, by = "marker") %>%
  select(trait, marker, chrom:cM_pos, trait, names(.)) %>%
  ungroup()

# Save the data
save("marker_mean_fw", file = file.path(result_dir, "marker_mean_fw_results.RData"))

