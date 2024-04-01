## Genomic prediction of FW sensitivity coefficients
## 
## Conduct genomic prediction and cross-validation of FW sensitivity coefficients
## For CV, use both GBS markers and BOPA markers. Also conduct predictions using the 
## genotypic means
## 
## Author: Jeff Neyhart
## Last updated: June 26, 2018
## 


# Run the source script - local
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

# Other packages
library(modelr)
library(parallel)
library(lme4)
library(broom)
library(cowplot)

# Load the FW results
load(file.path(result_dir, "pheno_mean_fw_results.RData"))
# Load the marker subsets
load(file.path(result_dir, "marker_subsets.RData"))

# Cores
n_cores <- detectCores()


# Rename the marker matrix
M <- S2TP_imputed_multi_genos_mat[tp_geno,]
# Overall K
K <- A.mat(X = M, min.MAF = 0, max.missing = 1)

# Colors for the coefficients
colors_use <- set_names(umn_palette(2)[3:5], coef_replace)


# Number of cv iterations
n_cv_iter <- 100


## TP cross-validation
p_train <- 0.60

# Create a tidy data.frame for modeling
# pheno_mean_fw_tomodel <- pheno_mean_fw %>%
pheno_mean_fw_tomodel <- pheno_mean_fw_tpvp %>%
  filter(line_name %in% tp) %>% # Do this to filter out only the TP
  distinct(trait, line_name, g, b, delta) %>%
  mutate(log_delta = log(delta)) %>%
  select(-delta) %>%
  gather(coef, value, g:log_delta)


## Create the cross-validation samples
pheno_mean_fw_tomodel_cv <- pheno_mean_fw_tomodel %>% 
  mutate(line_name = as.factor(line_name)) %>%
  group_by(trait, coef) %>% 
  do(cv_sample = crossv_mc(data = ., n = n_cv_iter, test = 1 - p_train)) %>%
  ungroup()

## Run cross-val
cv_results <- pheno_mean_fw_tomodel_cv %>%
  group_by(trait, coef) %>%
  do({
    
    df <- .
    
    # Map, predict, and measure accuracy
    acc_out <- df$cv_sample[[1]] %>% 
      pmap(~predict_RR(train = ..1, test = ..2, K = K)) %>%
      map_dbl(~cor(.$value, .$pred_value))
    
    df %>% 
      unnest() %>%
      mutate(iter = seq(n()), acc = acc_out) %>% 
      select(-train:-.id)
    
  })
  

# Compare correlations
cv_results %>% 
  summarize_at(vars(acc), funs(mean = mean(.), sd = sd(.)))


# Boxplot
g_cv_results_box <- cv_results %>% 
  ungroup() %>%
  mutate(coef = str_replace_all(coef, coef_replace)) %>%
  ggplot(aes(x = trait, y = acc, fill = coef)) + 
  geom_boxplot(position = "dodge") + 
  scale_fill_manual(values = colors_use, name = NULL) +
  ylab("Prediction accuracy") +
  ylim(c(-0.3, 1)) +
  theme_pnas() + 
  theme(axis.title.x = element_blank(), legend.position = c(0.80, 0.85))
  
# Save
ggsave(filename = "cv_all_boxplot.jpg", plot = g_cv_results_box, path = fig_dir,
       height = 9, width = 8, units = "cm", dpi = 1000)

    



## Run similar 60:40 CV, but instead predict the value of the unobserved genotypes in
## each environment

## This will assume two things: 1) the estimate of the mean within an environment is not 
## dependent on the population observed in that environment (reasonable assumption), and 2) that
## the estimates of linear stability do not depend on all environments (also reasonable).
pheno_mean_fw_tomodel <- pheno_mean_fw_tpvp %>%
  filter(line_name %in% tp) %>%
  select(trait, line_name, environment, value, g, b, h)


## Create the cross-validation samples
pheno_mean_fw_tomodel_cv <- pheno_mean_fw_tomodel %>% 
  mutate(line_name = as.factor(line_name)) %>%
  group_by(trait) %>% 
  do(cv_sample = {
    df <- . 
    
    # Create a set of training lines (we assume we know the mean and regression slope)
    sample_train_line_name <- replicate(n = n_cv_iter, sample_frac(tbl = distinct(df, line_name), size = p_train), simplify = FALSE)
    
    # Create training sets of the mean and stability of the TP
    training_data <- sample_train_line_name %>%
      map(~left_join(., distinct(df, trait, line_name, g, b), by = c("line_name")))
    
    # Create testing sets of the environment-specific observations of the testing pop
    testing_data <- sample_train_line_name %>%
      map(~anti_join(df, ., by = c("line_name")))
    
    # Return the training and testing sets
    data_frame(train = training_data, test = testing_data)
    
  })
    

## Run cross-val
cv_results <- pheno_mean_fw_tomodel_cv %>%
  group_by(trait) %>%
  do({
    df <- .
    
    # Iterate over cv iterations
    acc_out_list <- df$cv_sample[[1]] %>%
      pmap(~{
        # .x is the training set. use it to calculate marker effects for both the mean and the slope
        g_pred_out <- mixed.solve(y = .x$g, model.matrix(~ -1 + line_name, .x), K = K)$u %>%
          data_frame(line_name = names(.), pred_g = .)
        b_pred_out <- mixed.solve(y = .x$b, model.matrix(~ -1 + line_name, .x), K = K)$u %>%
          data_frame(line_name = names(.), pred_b = .)
        
        # Predict the value of the testing set in each environment based on the environmental mean
        list(.y, g_pred_out, b_pred_out) %>% 
          reduce(left_join, by = "line_name") %>% 
          mutate(pred_value_gb = pred_g + (pred_b * h),
                 pred_value_g = pred_g) %>% # Also just use the genotype mean 
          group_by(environment) %>% 
          summarize(acc_gb = cor(value, pred_value_gb),
                    acc_g = cor(value, pred_value_g))
        
      })
    
    # Add the accuracy results to the df
    df %>% 
      unnest() %>% 
      mutate(acc_out = acc_out_list, iter = seq(n())) %>% 
      select(trait, iter, acc_out)
    
  })
   
 
# Look at environment-specific distributions in prediction accuracy
cv_results_env_spec <- cv_results %>% 
  ungroup() %>% 
  unnest()

# Plot
g_cv_env_obs <- cv_results_env_spec %>%
  gather(pred_type, acc, acc_gb, acc_g) %>%
  mutate(pred_type = str_replace_all(pred_type, c(acc_gb = "MeanEnvPrediction", acc_g = "MeanPrediction"))) %>%
  ggplot(aes(x = environment, y = acc, fill = pred_type)) + 
  geom_boxplot() +
  facet_grid(trait ~ ., switch = "y", scales = "free_y") +
  scale_fill_manual(values = unname(colors_use)[1:2], name = NULL) +
  ylab("Prediction accuracy") +
  theme_pnas() +
  theme(axis.text.x = element_text(angle = 90), axis.title.x = element_blank(),
        panel.spacing.y = unit(1, "lines"), legend.position = "bottom")

## Save
ggsave(filename = "cv_env_spec_obs_boxplot.jpg", plot = g_cv_env_obs, path = fig_dir,
       height = 12, width = 17.4, units = "cm", dpi = 1000)



## There is a clear environment effect.
# Fit a model to get the mean accuracy
cv_results_env_spec_fit <- cv_results_env_spec %>%
  mutate(iter = as.factor(iter)) %>%
  group_by(trait) %>%
  do({
    df <- .
    fit <- lm(acc_gb ~ 1 + iter + environment, data = df, contrasts = list(environment = "contr.sum", iter = "contr.sum"))
    # Get the effects of each iteration (add the intercept to get the mean accuracy in each iteration)
    iter_coef <- coef(fit)[2:n_cv_iter]
    iter_effects <- c(iter_coef, c(iter100 = -sum(iter_coef)))
    iter_means <- iter_effects + coef(fit)[1]
    
    # Return a data.frame
    data_frame(iter = names(iter_means), acc = iter_means) %>%
      mutate(iter = parse_number(iter))
    
  })
    

## Simply calculate environment means over iterations
cv_results_env_summ <- cv_results_env_spec %>%
  group_by(trait, environment) %>%
  summarize_at(vars(acc_gb, acc_g), mean)
  



# Summarize the accuracy accross all environments for each iteration
cv_results_summ <- cv_results_env_spec  %>% 
  group_by(trait, iter) %>% 
  summarize_at(vars(acc_gb, acc_g), mean)

# Boxplot
g_cv_results_env_spec_box <- cv_results_summ %>% 
  ungroup() %>%
  gather(pred_type, acc, acc_gb:acc_g) %>%
  mutate(pred_type = str_replace_all(pred_type, c(acc_gb = "MeanEnvPrediction", acc_g = "MeanPrediction"))) %>%
  ggplot(aes(x = trait, y = acc, fill = pred_type)) + 
  geom_boxplot(position = "dodge") + 
  scale_fill_manual(values = unname(colors_use)[1:2], name = NULL) +
  ylab("Prediction accuracy") +
  ylim(c(0, 1)) +
  theme_pnas() + 
  theme(axis.title.x = element_blank(), legend.position = c(0.80, 0.85))

# Save
ggsave(filename = "cv_env_spec_boxplot.jpg", plot = g_cv_results_env_spec_box, path = fig_dir,
       height = 9, width = 8, units = "cm", dpi = 1000)



## Compare the distributions of these means using a Mann-Whitney test
# The test is paired because it is the same TP that is used for prediction in each iteration
cv_results_summ %>% 
  do(t_test = t.test(x = .$acc_gb, y = .$acc_g, paired = TRUE, conf.int = TRUE)) %>%
  ungroup() %>%
  mutate(t_test_pvalue = map_dbl(t_test, "p.value"))





#### Predictions of the VP


# Rename the marker matrix
M <- s2_imputed_mat[c(tp_geno1, vp_geno1),]
# Overall K
K <- A.mat(X = M, min.MAF = 0, max.missing = 1)


# Create a tidy dataset to model
pheno_mean_fw_tomodel <- pheno_mean_fw_tpvp %>% 
  filter(line_name %in% c(tp_geno1, vp_geno1)) %>%
  distinct(trait, line_name, g, b, delta) %>%
  mutate(line_name = as.factor(line_name),
         log_delta = log(delta)) %>%
  select(-delta) %>%
  gather(coef, value, g:log_delta)


## For each trait, predict the mean, slope, and MSE of the validation set
vp_prediction_all_data <- pheno_mean_fw_tomodel %>% 
  mutate(line_name = as.factor(line_name)) %>%
  group_by(trait, coef) %>%
  do({
    df <- .
    
    # Separate the training from validation
    df_train <- filter(df, line_name %in% tp_geno1)
    mf <- model.frame(value ~ line_name, df_train)
    y <- model.response(mf)
    Z <- model.matrix(~ -1 + line_name, mf)
    
    # Fit the model
    fit <- mixed.solve(y = y, Z = Z, K = K)
    
    # Return the BLUPs
    blups <- data.frame(line_name = names(fit$u), pred_value = fit$u, row.names = NULL, stringsAsFactors = FALSE)
    # Just the vp
    filter(df, line_name %in% vp_geno1) %>%
      left_join(., blups, by = "line_name")
    
  })


# Correlate
vp_prediction_summ <- vp_prediction_all_data %>% 
  summarize(acc = cor(value, pred_value)) %>%
  mutate(annotation = str_c("r[MG]~'= ", formatC(acc, digits = 3, format = "f"), "'"),
         coef = str_replace_all(coef, coef_replace))

spread(vp_prediction_summ, coef, acc)

# trait            b     g log_delta
# 1 GrainYield  0.289  0.553    0.0494
# 2 HeadingDate 0.472  0.535    0.215 
# 3 PlantHeight 0.0944 0.614    0.297 


# Plot
g_all_data_predict <- vp_prediction_all_data %>% 
  ungroup() %>%
  mutate(coef = str_replace_all(coef, coef_replace)) %>%
  ggplot(aes(x = value, pred_value, color = coef, shape = coef)) + 
  geom_point(size = 0.5) + 
  geom_smooth(method = "lm", se = FALSE, lwd = 0.5) + 
  geom_text(data = mutate(vp_prediction_summ, coef = str_replace_all(coef, coef_replace)), parse = TRUE,
            aes(x = Inf, y = -Inf, label = annotation), inherit.aes = FALSE, size = 1.5, hjust = 1, vjust = -1) +
  facet_wrap(trait ~ coef, scales = "free") + 
  scale_color_manual(values = colors_use, guide = FALSE) + 
  scale_shape_discrete(guide = FALSE) + 
  ylab("Predicted value") +
  xlab("Observed value") +
  theme_pnas()

# Save the plot
ggsave(filename = "vp_prediction_all_data.jpg", plot = g_all_data_predict, path = fig_dir,
       height = 9, width = 8.7, units = "cm", dpi = 1000)





## Create a plot per coef, then combine
# Prepare a df
data_toplot <- vp_prediction_all_data %>%
  ungroup() %>%
  mutate(coef = str_replace_all(coef, coef_replace),
         coef = as.factor(coef))

# Mean
g_all_data_predict_g <- data_toplot %>% 
  filter(coef == "Genotype Mean") %>%
  ggplot(aes(x = value, pred_value, color = coef, shape = coef)) + 
  geom_point(size = 0.5) + 
  geom_smooth(method = "lm", se = FALSE, lwd = 0.5) + 
  geom_text(data = filter(vp_prediction_summ, coef == "Genotype Mean"), parse = TRUE,
            aes(x = Inf, y = -Inf, label = annotation), inherit.aes = FALSE, size = 1.5, hjust = 1, vjust = -1) +
  facet_wrap(~ trait, scales = "free", strip.position = "left", ncol = 1) + 
  scale_color_manual(values = colors_use, guide = FALSE) + 
  scale_shape_discrete(guide = FALSE) + 
  ylab("Predicted value") +
  xlab("Observed value") +
  theme_pnas()

# Linear stability
g_all_data_predict_b <- data_toplot %>% 
  filter(coef == "Linear Stability") %>%
  ggplot(aes(x = value, pred_value, color = coef, shape = coef)) + 
  geom_point(size = 0.5) + 
  geom_smooth(method = "lm", se = FALSE, lwd = 0.5) + 
  geom_text(data = filter(vp_prediction_summ, coef == "Linear Stability"), parse = TRUE,
            aes(x = Inf, y = -Inf, label = annotation), inherit.aes = FALSE, size = 1.5, hjust = 1, vjust = -1) +
  facet_wrap(~ trait, scales = "free", strip.position = "left", ncol = 1) + 
  scale_color_manual(values = colors_use, name = NULL, drop = FALSE) + 
  scale_shape_discrete(guide = FALSE) + 
  ylab("Predicted value") +
  xlab("Observed value") +
  theme_pnas() +
  theme(legend.position = "bottom")

# Non-linear stability
g_all_data_predict_d <- data_toplot %>% 
  filter(coef == "Non-Linear Stability") %>%
  ggplot(aes(x = value, pred_value, color = coef, shape = coef)) + 
  geom_point(size = 0.5) + 
  geom_smooth(method = "lm", se = FALSE, lwd = 0.5) + 
  geom_text(data = filter(vp_prediction_summ, coef == "Non-Linear Stability"), parse = TRUE,
            aes(x = Inf, y = -Inf, label = annotation), inherit.aes = FALSE, size = 1.5, hjust = 1, vjust = -1) +
  facet_wrap(~ trait, scales = "free", strip.position = "left", ncol = 1) + 
  scale_color_manual(values = colors_use, guide = FALSE) + 
  scale_shape_discrete(guide = FALSE) + 
  ylab("Predicted value") +
  xlab("Observed value") +
  theme_pnas()


## Combine plots
g_vp_all_data_predict <- plot_grid(
  g_all_data_predict_g + theme(axis.title.x = element_blank(), strip.placement = "outside"),
  g_all_data_predict_b + theme(axis.title.y = element_blank(), strip.background = element_blank(),
                               strip.text = element_blank(), legend.position = "none"),
  g_all_data_predict_d + theme(axis.title = element_blank(), strip.background = element_blank(),
                               strip.text = element_blank()),
  ncol = 3, align = "hv", axis = "tb", rel_widths = c(1, 0.8, 0.8)
)

# Add legend
g_vp_all_data_predict1 <- plot_grid(
  g_vp_all_data_predict,
  get_legend(g_all_data_predict_b),
  ncol = 1, rel_heights = c(1, 0.15)
)

# Save
ggsave(filename = "vp_prediction_all_data_merge.jpg", plot = g_vp_all_data_predict1, path = fig_dir,
       height = 7, width = 8.7, units = "cm", dpi = 1000)










## Conduct the same predictions of environment-specific genotypes
## Again, this assumes that the environmental mean is known

# Create a tidy dataset to model
pheno_mean_fw_tomodel <- pheno_mean_fw_tpvp %>% 
  filter(line_name %in% c(tp_geno1, vp_geno1)) %>%
  mutate(line_name = as.factor(line_name)) %>%
  select(trait, line_name, environment, value, g, h, b)

cv_results <- pheno_mean_fw_tomodel %>%
  group_by(trait) %>%
  do({
    df <- .
    
    # Extract the training set
    train <- distinct(df, line_name, g, b) %>%
      filter(line_name %in% tp_geno1)
    
    # Predict the mean and the stability
    g_pred_out <- mixed.solve(y = train$g, model.matrix(~ -1 + line_name, train), K = K)$u %>%
      data_frame(line_name = names(.), pred_g = .)
    b_pred_out <- mixed.solve(y = train$b, model.matrix(~ -1 + line_name, train), K = K)$u %>%
      data_frame(line_name = names(.), pred_b = .)
    
    # Extract the test set
    test <- df %>% 
      filter(line_name %in% vp_geno1)
      
    list(test, g_pred_out, b_pred_out) %>% 
      reduce(left_join, by = "line_name") %>% 
      mutate(pred_value_gb = pred_g + (pred_b * h),
             pred_value_g = pred_g) %>% # Also just use the genotype mean 
      group_by(environment) %>% 
      summarize(acc_gb = cor(value, pred_value_gb),
                acc_g = cor(value, pred_value_g))
    
  })


# Summarize the accuracy accross all environments for each iteration
cv_results_summ <- cv_results  %>% 
  group_by(trait) %>% 
  summarize_at(vars(acc_gb, acc_g), mean)

# trait       acc_gb acc_g
# 1 GrainYield   0.243 0.232
# 2 HeadingDate  0.419 0.414
# 3 PlantHeight  0.411 0.410

# Boxplot
g_cv_results_env_spec_box <- cv_results %>% 
  ungroup() %>%
  gather(pred_type, acc, acc_gb:acc_g) %>%
  mutate(pred_type = str_replace_all(pred_type, c(acc_gb = "MeanEnvPrediction", acc_g = "MeanPrediction"))) %>%
  ggplot(aes(x = trait, y = acc, fill = pred_type)) + 
  geom_boxplot(position = "dodge") + 
  scale_fill_manual(values = unname(colors_use)[1:2], name = NULL) +
  ylab("Prediction accuracy") +
  ylim(c(0, 1)) +
  theme_pnas() + 
  theme(axis.title.x = element_blank(), legend.position = c(0.80, 0.85))


## Compare with a t-test
cv_results %>% 
  do(t_test = t.test(x = .$acc_gb, y = .$acc_g, paired = TRUE)) %>% 
  ungroup() %>% 
  mutate(t_test_pvalue = map_dbl(t_test, "p.value"))

# trait       t_test      t_test_pvalue
# 1 GrainYield  <S3: htest>        0.0828
# 2 HeadingDate <S3: htest>        0.113 
# 3 PlantHeight <S3: htest>        0.887








# VP years
vp_years <- unique(S2_MET_BLUEs_tpvp$year)

# Empty list
vp_prediction_loyo <- list()

# Iterate over the years
for (yr in vp_years) {
  
  # Extract the TP and VP data
  pheno_use_tp <- pheno_mean_fw_tpvp %>% 
    filter(line_name %in% tp_geno1, year != yr)
  
  pheno_use_vp <- pheno_mean_fw_tpvp %>% 
    filter(line_name %in% vp_geno1, year == yr)
  
  ## Calculate g and h
  tp_gh <- pheno_use_tp %>%
    select(trait, line_name, environment, value, std_error) %>%
    group_by(trait) %>%
    do(calc_gh(.))
  
  # Calculate stability
  tp_stab <- tp_gh %>%
    group_by(trait, line_name) %>%
    do(calc_stability(., remove.outliers = FALSE)) %>%
    filter(type == "outliers") %>%
    select(-type:-model)
  
  # Combine
  tp_pheno_fw <- left_join(tp_gh, tp_stab)
  
  ## Do the same for the vp
  vp_gh <- pheno_use_vp %>%
    select(trait, line_name, environment, value, std_error) %>%
    group_by(trait) %>%
    do(calc_gh(.))
  
  # Calculate stability
  vp_stab <- vp_gh %>%
    group_by(trait, line_name) %>%
    do(calc_stability(., remove.outliers = FALSE)) %>%
    filter(type == "outliers") %>%
    select(-type:-model)
  
  # Combine
  vp_pheno_fw <- left_join(vp_gh, vp_stab)
  
  
  
  ## Predict
  # Create a data.frame with all data
  pheno_tomodel <- bind_rows(tp_pheno_fw, vp_pheno_fw) %>%
    distinct(trait, line_name, g, b, delta) %>%
    mutate(log_delta = log(delta)) %>%
    select(-delta) %>%
    gather(coef, value, g:log_delta)
  
  vp_prediction <- pheno_tomodel %>% 
    mutate(line_name = as.factor(line_name)) %>%
    group_by(trait, coef) %>%
    do({
      df <- .
      
      # Separate the training from validation
      df_train <- filter(df, line_name %in% tp_geno1)
      mf <- model.frame(value ~ line_name, df_train)
      y <- model.response(mf)
      Z <- model.matrix(~ -1 + line_name, mf)
      
      # Fit the model
      fit <- mixed.solve(y = y, Z = Z, K = K)
      
      # Return the BLUPs
      blups <- data.frame(line_name = names(fit$u), pred_value = fit$u, row.names = NULL, stringsAsFactors = FALSE)
      # Just the vp
      filter(df, line_name %in% vp_geno1) %>%
        mutate(line_name = as.character(line_name)) %>%
        left_join(., blups, by = "line_name")
      
    }) %>%
    mutate(n_vp_env = n_distinct(vp_pheno_fw$environment),
           year = yr)
  
  
  # Add to the list
  vp_prediction_loyo[[as.character(yr)]] <- vp_prediction
  
}


## Consolidate
vp_prediction_loyo_acc <- vp_prediction_loyo %>%
  map(~summarize(., acc = cor(value, pred_value))) %>% 
  list(., names(.)) %>% 
  pmap_df(~mutate(.x, vp_year = .y))






## Randomly subset environments 75%-25%
p_env_train <- 0.70
n_iter <- 25

# Generate environmental subsets
envs <- unique(S2_MET_BLUEs_tpvp$environment)
tp_env_samples <- replicate(n = n_iter, sample(x = envs, size = floor(p_env_train * length(envs))), simplify = FALSE)

# Iterate over the samples
vp_prediction_rand_env <- tp_env_samples %>%
  map(function(sample_envs) {
    
    # Extract the TP and VP data
    pheno_use_tp <- pheno_mean_fw_tpvp %>% 
      filter(line_name %in% tp_geno1, environment %in% sample_envs)
    
    pheno_use_vp <- pheno_mean_fw_tpvp %>% 
      filter(line_name %in% vp_geno1, !environment %in% sample_envs)
    
    ## Calculate g and h
    tp_gh <- pheno_use_tp %>%
      select(trait, line_name, environment, value, std_error) %>%
      group_by(trait) %>%
      do(calc_gh(.))
    
    # Calculate stability
    tp_stab <- tp_gh %>%
      group_by(trait, line_name) %>%
      do(calc_stability(., remove.outliers = FALSE)) %>%
      filter(type == "outliers") %>%
      select(-type:-model)
    
    # Combine
    tp_pheno_fw <- left_join(tp_gh, tp_stab)
    
    ## Do the same for the vp
    vp_gh <- pheno_use_vp %>%
      select(trait, line_name, environment, value, std_error) %>%
      group_by(trait) %>%
      do(calc_gh(.))
    
    # Calculate stability
    vp_stab <- vp_gh %>%
      group_by(trait, line_name) %>%
      do(calc_stability(., remove.outliers = FALSE)) %>%
      filter(type == "outliers") %>%
      select(-type:-model)
    
    # Combine
    vp_pheno_fw <- left_join(vp_gh, vp_stab)
    
    
    
    ## Predict
    # Create a data.frame with all data
    pheno_tomodel <- bind_rows(tp_pheno_fw, vp_pheno_fw) %>%
      distinct(trait, line_name, g, b, delta) %>%
      mutate(log_delta = log(delta)) %>%
      select(-delta) %>%
      gather(coef, value, g:log_delta)
    
    vp_prediction <- pheno_tomodel %>% 
      mutate(line_name = as.factor(line_name)) %>%
      group_by(trait, coef) %>%
      do({
        df <- .
        
        # Separate the training from validation
        df_train <- filter(df, line_name %in% tp_geno1)
        mf <- model.frame(value ~ line_name, df_train)
        y <- model.response(mf)
        Z <- model.matrix(~ -1 + line_name, mf)
        
        # Fit the model
        fit <- mixed.solve(y = y, Z = Z, K = K)
        
        # Return the BLUPs
        blups <- data.frame(line_name = names(fit$u), pred_value = fit$u, row.names = NULL, stringsAsFactors = FALSE)
        # Just the vp
        filter(df, line_name %in% vp_geno1) %>%
          mutate(line_name = as.character(line_name)) %>%
          left_join(., blups, by = "line_name")
        
      }) %>%
      mutate(n_vp_env = n_distinct(vp_pheno_fw$environment))
    
  })


## Consolidate
vp_prediction_rand_env_acc <- vp_prediction_rand_env %>% 
  map(~summarize(., acc = cor(value, pred_value))) %>% 
  list(., seq_along(.)) %>% 
  pmap_df(~mutate(.x, iter = .y))





## Use marker subsets to predict stability  
# Load the marker subsets
load(file.path(result_dir, "marker_subsets_tpvp.RData"))

## Top markers
top_rank_markers_predictions <- top_rank_markers_tpvp %>%
  map(function(marker_df) {
    
    pheno_mean_fw_tomodel %>% 
      mutate(line_name = as.factor(line_name)) %>%
      group_by(trait, coef) %>%
      do({
        df <- .
        
        # Grab the markers
        markers_use <- subset(marker_df, trait == unique(df$trait) & coef == unique(df$coef), marker, drop = TRUE)
        # Create relationship matrix
        K_mat <- A.mat(X = M[,markers_use], min.MAF = 0, max.missing = 1)
        
        # Separate the training from validation
        df_train <- filter(df, line_name %in% tp_geno1)
        mf <- model.frame(value ~ line_name, df_train)
        y <- model.response(mf)
        Z <- model.matrix(~ -1 + line_name, mf)
        
        # Fit the model
        fit <- mixed.solve(y = y, Z = Z, K = K_mat)
        
        # Return the BLUPs
        blups <- data.frame(line_name = names(fit$u), pred_value = fit$u, row.names = NULL, stringsAsFactors = FALSE)
        # Just the vp
        filter(df, line_name %in% vp_geno1) %>%
          left_join(., blups, by = "line_name")
        
      })
    
  })
    


## Top, evenly-spaced markers
tesm_markers_predictions <- top_rank_evenly_spaced_markers_tpvp %>%
  map(function(marker_df) {
    
    pheno_mean_fw_tomodel %>% 
      mutate(line_name = as.factor(line_name)) %>%
      group_by(trait, coef) %>%
      do({
        df <- .
        
        # Grab the markers
        markers_use <- subset(marker_df, trait == unique(df$trait) & coef == unique(df$coef), marker, drop = TRUE)
        # Create relationship matrix
        K_mat <- A.mat(X = M[,markers_use], min.MAF = 0, max.missing = 1)
        
        # Separate the training from validation
        df_train <- filter(df, line_name %in% tp_geno1)
        mf <- model.frame(value ~ line_name, df_train)
        y <- model.response(mf)
        Z <- model.matrix(~ -1 + line_name, mf)
        
        # Fit the model
        fit <- mixed.solve(y = y, Z = Z, K = K_mat)
        
        # Return the BLUPs
        blups <- data.frame(line_name = names(fit$u), pred_value = fit$u, row.names = NULL, stringsAsFactors = FALSE)
        # Just the vp
        filter(df, line_name %in% vp_geno1) %>%
          left_join(., blups, by = "line_name")
        
      })
    
  })


## Evenly-spaced markers
esm_markers_predictions <- evenly_spaced_markers_tpvp %>%
  map(function(marker_df) {
    
    # Grab the markers
    markers_use <- subset(marker_df, , marker, drop = TRUE)
    # Create relationship matrix
    K_mat <- A.mat(X = M[,markers_use], min.MAF = 0, max.missing = 1)
    
    pheno_mean_fw_tomodel %>% 
      mutate(line_name = as.factor(line_name)) %>%
      group_by(trait, coef) %>%
      do({
        df <- .
        
        # Separate the training from validation
        df_train <- filter(df, line_name %in% tp_geno1)
        mf <- model.frame(value ~ line_name, df_train)
        y <- model.response(mf)
        Z <- model.matrix(~ -1 + line_name, mf)
        
        # Fit the model
        fit <- mixed.solve(y = y, Z = Z, K = K_mat)
        
        # Return the BLUPs
        blups <- data.frame(line_name = names(fit$u), pred_value = fit$u, row.names = NULL, stringsAsFactors = FALSE)
        # Just the vp
        filter(df, line_name %in% vp_geno1) %>%
          left_join(., blups, by = "line_name")
        
      })
    
  })


## Random markers
# Create a df to parallelize
random_markers_split <- data_frame(nmar = names(random_markers_tpvp), random_markers_tpvp) %>% 
  unnest() %>% 
  group_by(nmar) %>%
  mutate(iter = seq(n())) %>% 
  ungroup() %>%
  assign_cores(n_cores) %>%
  split(.$core)

# Parallelize
rand_marker_predction_out <- mclapply(X = random_markers_split, FUN = function(core_df) {
  
  # Map over the markers
  out <- core_df$random_markers_tpvp %>%
    map("marker") %>%
    map(~{
      
      K_mat <- A.mat(X = M[,.], min.MAF = 0, max.missing = 1)
      
      pheno_mean_fw_tomodel %>% 
        mutate(line_name = as.factor(line_name)) %>%
        group_by(trait, coef) %>%
        do({
          df <- .
          
          # Separate the training from validation
          df_train <- filter(df, line_name %in% tp_geno1)
          mf <- model.frame(value ~ line_name, df_train)
          y <- model.response(mf)
          Z <- model.matrix(~ -1 + line_name, mf)
          
          # Fit the model
          fit <- mixed.solve(y = y, Z = Z, K = K_mat)
          
          # Return the BLUPs
          blups <- data.frame(line_name = names(fit$u), pred_value = fit$u, row.names = NULL, stringsAsFactors = FALSE)
          # Just the vp
          filter(df, line_name %in% vp_geno1) %>%
            left_join(., blups, by = "line_name")
          
        })
      
    })
  
  core_df %>%
    mutate(results = out) %>%
    select(-random_markers_tpvp, -core)
  
}, mc.cores = n_cores)


rand_marker_predctions <- bind_rows(rand_marker_predction_out)


## Save everything
all_marker_prediction_results <- list(
  vp_prediction_all_data = vp_prediction_all_data,
  vp_prediction_loyo = vp_prediction_loyo,
  vp_prediction_rand_env = vp_prediction_rand_env
)

marker_subset_prediction_results <- list(
  top_rank_markers_predictions = top_rank_markers_predictions,
  tesm_markers_predictions = tesm_markers_predictions,
  esm_markers_predictions = esm_markers_predictions,
  rand_marker_predctions = rand_marker_predctions
)

save("all_marker_prediction_results", "marker_subset_prediction_results", file = file.path(result_dir, "vp_stability_prediction.RData"))



