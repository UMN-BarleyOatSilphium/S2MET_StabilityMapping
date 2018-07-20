## S2MET Mapping Functions
## 
## A script with useful functions for the S2MET mapping project
## 



# A function to assign cores to a data.frame
assign_cores <- function(df, n_core) {
  df$core <- sort(rep(seq(n_core), length.out = nrow(df)))
  return(df)
}



## Calculate marginal and conditional R-squared values from the output of `mixed.solve`
r_squared <- function(object, snps, X) {
  
  # Calculated predicted values using all fixed effects
  fixed_pred <- X %*% object$beta
  # Get the random effect variance components (and residual)
  ranef_varcomp <- object$Vu
  resid_varcomp <- object$Ve
  
  # Calculate the denominator of the R-squared equations
  denom <- c(var(fixed_pred) + ranef_varcomp + resid_varcomp)
  
  # Now calculate the R-squared explained by each snp
  fixed_pred_each <- t(c(object$beta[snps]) * t(X[,snps, drop = FALSE]))
  fixed_var_each <- apply(X = fixed_pred_each, MARGIN = 2, FUN = var)
  
  r_sqrd_snp <- fixed_var_each / denom
  
  # Calculate R-squared for all SNPs jointly
  fixed_pred_snps <- X[,snps,drop = FALSE] %*% object$beta[snps] 
  r_sqrd_snps <- var(fixed_pred_snps) / denom
  
  ## All fixed effects
  r_sqrd_fixed <- var(fixed_pred) / denom
  ## Fixed and random
  r_sqrd_eff <- (var(fixed_pred) + ranef_varcomp) / denom
  
  # Return a list
  list(
    fixed_r_squared = c(all_fixed = c(r_sqrd_fixed), fixed_rand = c(r_sqrd_eff), all_snps = r_sqrd_snps),
    snp_r_squared = c(r_sqrd_snp)
  )

}



## A function to fit a multi-locus mixed model and stepwise remove
## non-significant loci.
mlmm <- function(y, X, Z, K, snps) {
  
  # Fit the model
  fit <- mixed.solve(y = y, Z = Z, K = K, X = X, SE = TRUE)
  beta <- fit$beta[-1]
  se <- fit$beta.SE[-1]
  # Perform hypothesis test
  chisq <- beta^2 / se^2
  pvalue <- pchisq(q = chisq, df = 1, lower.tail = FALSE)
  qvalue <- p.adjust(p = pvalue, method = "BH")
  
  # Get r_squared values
  R2 <- r_squared(object = fit, snps = snps, X = X)
  
  # Return a list and a data.frame
  fixed_sig <- data.frame(term = names(beta), beta, se, pvalue, qvalue, row.names = NULL)
  list(summary = fixed_sig, r_squared = R2)
  
}


# Function for removing colinear SNPs
remove_colinear <- function(x) {
  # LD
  x_LD <- cor(x)^2
  # Find the index of the SNPs to be dropped
  snp_index <- apply(X = x_LD, MARGIN = 1, FUN = function(LD) max(which(LD > 0.99)))
  to_drop <- which(duplicated(snp_index))
  
  ifelse(length(to_drop) > 0, return(x[,-to_drop, drop = FALSE]), return(x))
  
}



## Multivariate association mapping for pleitropy
mv_gwas <- function(pheno, geno, fixed = NULL, K, n.PC) {
  
  # Stop if pheno does not have at least 3 cols
  stopifnot(ncol(pheno) >= 3)
  
  # Get the line names and the name of the column storing that information
  linnenames <- pheno[,1]
  linname_col <- colnames(pheno)[1]
  # Subset the SNP data
  snp_info <- geno[,1:3]
  
  # Create the genotype matrix and add the marker names
  M <- t(geno[,linnenames,drop = FALSE])
  colnames(M) <- snp_info[,1]
  
  # Calculate K, if absent
  if (missing(K)) K <- A.mat(X = M, min.MAF = 0, max.missing = 1)
  
  # If PC > 1, get the eigenvector
  if (n.PC > 0) {
    K_eigen <- eigen(K)$vectors
    PCs <- K_eigen[,seq(n.PC),drop = FALSE]
    colnames(PCs) <- paste0("PC", seq(n.PC))
    
  } else {
    PCs <- NULL
    
  }
  
  ## Fit the mutli-variate model
  # First create the model matrices
  Y <- as.matrix(pheno[,-1])
  # Number of parameters
  d <- ncol(Y)
  # Trait names
  traits <- colnames(Y)
  
  print(paste("Multivariate GWAS for trait pair:", paste(traits, collapse = ", ")))

  # Create the fixed effect matrix
  mu <- model.matrix(~ 1, pheno)
  X <- cbind(mu, PCs)
  
  # Random effect matrix
  form <- as.formula(paste("~ -1 +", linname_col))
  Z <- model.matrix(form, pheno)
  
  ## Fit the model
  fit_base <- emmremlMultivariate(Y = t(Y), X = t(X), Z = t(Z), K = K)
  # Calculate the Hinv matrix
  ZKZt <- (Z %*% K %*% t(Z))
  Hinv_mv <- solve(kronecker(ZKZt, fit_base$Vg) + kronecker(diag(nrow(ZKZt)), fit_base$Ve))
  
  # Vectorize the Y matrix
  Yvec <- Y[rep(seq(nrow(Y)), each = d),,drop = FALSE]
  mu_vec <- kronecker(mu, diag(d))
  
  # Number of samples and parameters
  n <- length(Y)
  p <- ncol(Y) * (ncol(X) + 1)
  v2 <- n - p
  seq_p <- setdiff(seq(p), seq(ncol(X) * d))
  
  
    
  print("Variance components estimated. Testing markers.")
  
  ## Iterate over the markers
  marker_score <- apply(X = M, MARGIN = 2, FUN = function(snp) {
    
    # Create a new X matrix
    X1 <- cbind(X, snp)
    # Replicate
    X1_use <- X1[rep(seq(nrow(X1)), each = d),,drop = FALSE]
    # Vectorize the X1
    Xforvec <- (kronecker(X1, diag(d)))
    
    W <- crossprod(X1_use, Hinv_mv %*% X1_use)
    Winv <- try(solve(W), silent = TRUE)
    if (class(Winv) != "try-error") {
      beta <- t(Winv %*% crossprod(X1_use, Hinv_mv %*% Yvec))
      CovBeta <- solve(crossprod(Xforvec, Hinv_mv %*% Xforvec))
      
      # Conduct Wald test
      statistic <- as.vector(beta[seq_p])^2 / diag(CovBeta)[seq_p]
      pvalue <- pchisq(q = statistic, df = 1, lower.tail = FALSE)
    } else {
      pvalue <- NA
    }
    
    return(pvalue) })
  
  # Transpose and add column names
  marker_score1 <- -log10(t(marker_score))
  colnames(marker_score1) <- colnames(Y)
  
  # Combine with the SNP information and return
  cbind(snp_info, marker_score1)

}


# Return an LD data.frame from a SNP matrix
LD <- function(x, df = TRUE) {
  ld <- cor(x)^2
  
  if (df) {
    ld <- as.data.frame(ld) %>% 
      rownames_to_column("marker1") %>% 
      gather(marker2, LD, -marker1)
  }
  return(ld)
}




## Finlay-Wilkinson Regression
# Use the mean of each genotype in each environment and the random effect of 
# each environment to calculate the random regression coefficient of genotype
# on environment


# Create a function to iteratively remove outliers based on studentized residuals
# Returns a df to be used to fit the final model
remove_outliers <- function(df, fit, cutoff = 3) {
  
  # Add residuals to the df
  df_resid <- df %>% 
    add_residuals(fit) %>%
    mutate(stand_resid = as.numeric(scale(resid, center = FALSE)))
  
  # If there are std residuals that are greater than the cutoff, remove them
  # and refit
  while (any(abs(df_resid$stand_resid) > cutoff)) {
    
    # Filter
    df_filter <- df_resid %>%
      filter(abs(stand_resid) <= cutoff)
    
    # Refit
    re_fit <- lm(value ~ h, data = df_filter)
    
    # Add the residuals and repeat
    df_resid <- df_filter %>% 
      add_residuals(re_fit) %>%
      mutate(stand_resid = as.numeric(scale(resid, center = FALSE)))
  }
  
  # Return the data.frame
  return(df_resid)
}

# Function to estimate the slope and MSE of a line
# Incorportate code to remove outliers
calc_stability <- function(df, remove.outliers = TRUE) {
  
  # # Set the control
  # control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore", check.nlev.gtr.1 = "ignore")
  # # Get the weights
  # wts <- df$std_error^2
  # fit <- lmer(value ~ (h|line_name), data = df, control = control, weights = wts)
  
  # Fit the model and return
  fit <- lm(value ~ h, data = df) 
  
  if (remove.outliers) {
  
    # Add residuals to the df
    # Filter for outliers
    df_filter <- remove_outliers(df = df, fit = fit, cutoff = 3)
    
    # Refit
    re_fit <- lm(value ~ h, data = df_filter) 
    
    df_re_fit <- df_filter %>%
      select(-contains("resid")) %>%
      mutate(b = coef(re_fit)[2],
             b_std_error = subset(tidy(re_fit), term == "h", std.error, drop = TRUE), # Regression coefficient
             delta = mean(resid(re_fit)^2))
  
  } else {
    df_re_fit <- data_frame(b = NA, delta = NA)
    re_fit <- NA
    
  }
  

  
  # Return a data.frame with each data.frame, slope and MSE estimates,
  # and number of observations
  df_fit <- df %>%
    mutate(b = coef(fit)[2],
           delta = mean(resid(fit)^2))
  

  
  # List of data.frame
  df_list <- list(df_fit, df_re_fit)
  
  # Return data
  df1 <- data_frame(type = c("outliers", "no_outliers"), 
                    data = df_list, 
                    model = list(fit, re_fit),
                    n = map_dbl(df_list, nrow), 
                    b = map_dbl(df_list, ~unique(.$b)), 
                    delta = map_dbl(df_list, ~unique(.$delta)))
  
  return(df1)
  
}


## Genotype mean
# Function to calculate the genotype mean and the environmental effect
calc_gh <- function(df) {
  
  # Set the control
  control <- lmerControl(check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
  # Extract the weights - square of the standard error of the mean
  wts <- df$std_error^2
  
  # First calculate the least-square means of each genotype
  # Fit environment and GxE as random effects
  fit1 <- lmer(value ~ -1 + line_name + (1|environment),
               data = df, control = control, weights = wts)
  
  # Extract the BLUEs
  geno_blues <- tidy(fit1) %>% 
    filter(group == "fixed") %>% 
    mutate(line_name = str_replace(term, "line_name", "")) %>% 
    select(line_name, g = estimate)
  
  # Extract the BLUPs
  env_blups <- ranef(fit1) %>% 
    as.data.frame() %>% 
    filter(grpvar == "environment") %>% 
    select(environment = grp, h = condval)
  
  # Return
  df1 <- df %>% 
    left_join(., geno_blues, by = "line_name") %>% 
    left_join(., env_blups, by = "environment")
  
  return(df1)
  
}



# Write a function that takes a train and test set and predicts using rrBLUP	+
predict_RR <- function(train, test, K) {	
  
  # Convert to df	
  train_df <- as.data.frame(train)	
  test_df <- as.data.frame(test)	
  
  # Create the model matrix	
  mf <- model.frame(value ~ line_name, train_df)	
  y <- model.response(mf)	
  Z <- model.matrix(~ -1 + line_name, mf)	
  
  fit <- mixed.solve(y = y, Z = Z, K = K)	
  
  # Tidy	
  u_hat_tidy <- fit$u %>% 	
    data.frame(line_name = names(.), pred_value = ., stringsAsFactors = FALSE, row.names = NULL)	
  
  # Combine and return the predictions	
  suppressWarnings(left_join(test_df, u_hat_tidy, by = "line_name"))	
  
}



# Write a function that takes a train and test set and predicts using rrBLUP
predict_RR <- function(train, test, K) {
  
  # Convert to df
  train_df <- as.data.frame(train)
  test_df <- as.data.frame(test)
  
  # Create the model matrix
  mf <- model.frame(value ~ line_name, train_df)
  y <- model.response(mf)
  Z <- model.matrix(~ -1 + line_name, mf)
  
  fit <- mixed.solve(y = y, Z = Z, K = K)
  
  # Tidy
  u_hat_tidy <- fit$u %>% 
    data.frame(line_name = names(.), pred_value = ., stringsAsFactors = FALSE, row.names = NULL)
  
  # Combine and return the predictions
  suppressWarnings(left_join(test_df, u_hat_tidy, by = "line_name"))
  
}


