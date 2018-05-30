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


