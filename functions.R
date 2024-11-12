## S2MET Mapping Functions
## 
## A script with useful functions for the S2MET mapping project
## 



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
  fixed_pred_each <- t(c(object$beta[-1]) * t(X[,-1, drop = FALSE]))
  fixed_var_each <- apply(X = fixed_pred_each, MARGIN = 2, FUN = var)
  
  r_sqrd_fixed_each <- fixed_var_each / denom
  
  # Calculate R-squared for all SNPs jointly
  fixed_pred_snps <- X[,snps,drop = FALSE] %*% object$beta[snps] 
  r_sqrd_snps <- var(fixed_pred_snps) / denom
  
  ## All fixed effects
  r_sqrd_fixed <- var(fixed_pred) / denom
  ## Random effects
  r_sqrd_rand <- ranef_varcomp / denom
  ## Fixed and random
  r_sqrd_eff <- (var(fixed_pred) + ranef_varcomp) / denom
  
  # Return a list
  list(
    r_squared = c(all_fixed = c(r_sqrd_fixed), all_rand = c(r_sqrd_rand), fixed_and_rand = c(r_sqrd_eff), all_snps = r_sqrd_snps),
    fixed_r_squared = c(r_sqrd_fixed_each)
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
  n <- nrow(Y)
  
  # Standardize the kinship matrix
  K_stand <- (n-1)/ sum((diag(n) - matrix(1,n,n)/n) * K) * K
  
  
  print(paste("Multivariate GWAS for trait pair:", paste(traits, collapse = ", ")))

  # Create the fixed effect matrix
  mu <- model.matrix(~ 1, pheno)
  X <- cbind(mu, PCs)
  
  # Random effect matrix
  form <- as.formula(paste("~ -1 +", linname_col))
  Z <- model.matrix(form, pheno)
  colnames(Z) <- colnames(K)
  
  Y_use <- Y
  # Y_use <- scale(Y)
  
  K_use <- K
  # K_use <- K_stand
  
  # ## Fit the model
  # fit_base <- emmremlMultivariate(Y = t(Y_use), X = t(X), Z = t(Z), K = K_use)
  # varG <- fit_base$Vg
  # varR <- fit_base$Ve
  
  fit_base <- sommer::mmer(Y = Y_use, X = X, Z = list(gen = list(Z = Z, K = K_use)), check.model = FALSE, silent = TRUE, draw = FALSE)
  varG <- fit_base$var.comp$gen
  varR <- fit_base$var.comp$units
  
  # Extract variance components to estimate a scalar (for P3D)
  K_comb <- kronecker(varG, K_use)
  I_comb <- kronecker(varR, diag(nrow(Y_use)))
  bigK <- K_comb + I_comb
  Mmat <- solve(bigK)
  
  # Vectorize Y matrix
  Yvec <- as.vector(Y_use)
  # Incidence for mean
  X1 <- cbind(kronecker(diag(d), mu), kronecker(diag(d), PCs))
  
  # Number of samples and parameters
  n <- nrow(Y)
  p <- ncol(Y) + (ncol(Y) * ncol(PCs)) + 1 # Number of intercepts + beta_PC + beta_snp
  v2 <- n - p

  marker_score <- apply(X = M, MARGIN = 2, FUN = function(snp) {
  
    ## Combine with SNPs
    X2 <- cbind(X1, snp = rep(snp, d))
    
    W <- crossprod(X2, Mmat %*% X2)
    Winv <- try(solve(W), silent = TRUE)
    
    if (class(Winv) != "try-error") {
      
      # Stolen from GWAS function in rrBLUP
      beta <- t(Winv %*% crossprod(X2, Mmat %*% Yvec))
      resid <- Yvec - X2 %*% t(beta)
      
      s2 <- as.double(crossprod(resid, Mmat %*% resid))/v2
      CovBeta <- s2 * Winv
      
      # Wald test
      Fstat <- beta[p]^2/CovBeta[p, p]
      pvalue <- pchisq(q = Fstat, df = 1, lower.tail = FALSE)
    } else {
      pvalue <- NA
    }
    
    return(pvalue)
    
  })
  
  
  # Transpose and add column names
  marker_score1 <- -log10(marker_score)
  
  # Combine with the SNP information and return
  out <- cbind(snp_info, marker_score1)
  names(out)[ncol(out)] <- paste0(colnames(Y), collapse = ".")
  
  return(out)
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




## For a set of pvalues and qvalues, return an FDR-corrected threshold
pvalue_fdr <- function(p, q, fdr = 0.05, window = 30) {
  
  # ## For testing
  # p <- 10^-df$neg_log_p
  # q <- df$qvalue
  
  # Combine p and q values
  values <- cbind(p, q)
  values1 <- values[order(q),]
  
  # Is the lowest q-values less than the threshold?
  if (values1[1,2] <= fdr) {
    
    # If so, find the q value that is closest to the fdr level
    x <- which.min(abs(values1[,2] - fdr))
    
    # Create a window around this value
    wind <- window
    
    first <- max(1, x - ceiling(wind / 2))
    last <- min(x + ceiling(wind / 2), nrow(values1))
    if ((last - first) < wind) {
      last <- last + (wind - last - first)
    }

    # Reverse the columns
    values2 <- values1[first:last,c(2,1)]
    # values2 <- values1[,c(2,1)]

    # Fit a spline curve and predict where the FDR level would fall
    spl <- smooth.spline(x = values2, df = 3)
    
    ## Find the low point
    # Range in q values
    rq <- range(values2[,1])
    pred_p <- predict(spl, seq(rq[1], rq[2], by = 0.001))$y
    
    pfdr <- predict(spl, x = fdr)$y
    
    ## If this is negative, subtract the minimum fitted value
    ifelse(pfdr < 0, pfdr - min(pred_p), pfdr)

  } else {
    # Otherwise return NA
    return(NA)
  }
  
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


## Genotype mean
# Function to calculate the genotype mean and the environmental effect
calc_gh <- function(object) {
  
  # Extract environmental effects
  fixefs <- fixef(object)
  
  # Get the levels of environment
  mf <- model.frame(object)
  env_levels <- levels(mf$environment)
  
  # Environment effects
  env_effects <- tail(fixefs, -1)
  env_effects <- c(env_effects, -sum(env_effects))
  env_effects <- unname(env_effects)
  
  # Return a tibble
  tibble(environment = env_levels, h = env_effects)
  
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

gblup_crossv <- function(fixed, random, K, data = NULL, .train = NULL, .test = NULL, method = c("rrBLUP", "sommer")) {
  # Error
  stopifnot(inherits(fixed, "formula"))
  stopifnot(inherits(random, "formula"))
  if (!is.null(data)) stopifnot(is.data.frame(data))
  if (!is.null(.train)) {
    .train <- as.data.frame(.train)
  }
  if (!is.null(.test)) {
    .test <- as.data.frame(.test)
  }
  method <- match.arg(method)
  
  # Try to find the random term in .train
  text <- attr(terms(random), "term.labels")
  regex_random <- regexpr(pattern = paste0(names(.train), collapse = "|"), text = text)
  random_term <- substring(text = text, first = regex_random, regex_random + attr(regex_random, "match.length") - 1)
  
  
  ### 
  ### Save this section for data
  ### 
  
  # Fit using .train, validate using .test
  if (method == "rrBLUP") {
    # Check if the random formula as 1 term
    if (length(attr(terms(random), "term.labels")) > 1) {
      stop("Only 1 random effect term is allowed when 'method' = 'rrBLUP'.")
    }
    
    # Check if fixed is multivariate
    if (ncol(as.matrix(model.response(model.frame(formula = fixed, data = .train)))) > 1) {
      stop("Multivariate models are not permitted when 'method' = 'rrBLUP'.")
    }
    

    fixed_term <- attr(terms(fixed), "term.labels")
    fixed_term <- if (length(fixed_term) == 0) NULL else fixed_term
    
    # Fit using kin.blup
    fit <- kin.blup(data = as.data.frame(.train), geno = random_term, pheno = as.character(fixed)[2], 
                    K = K, fixed = fixed_term)
    
    # Get the predictions
    preds <- tibble(term = names(fit$pred), predicted.value = as.numeric(fit$pred))
    names(preds)[1] <- random_term
    
  } else if (method == "sommer") {
    # Fit with sommer
    fit <- mmer(fixed = fixed, random = random, data = .train, verbose = FALSE)
    ranef <- randef(object = fit)
    fixef <- as.matrix(fit$Beta[,"Estimate", drop = FALSE])
    
    # respone names
    response_names <- fit$terms$response[[1]]
    # Get the BLUPS
    blup_mat <- as.matrix(as.data.frame(ranef[[paste0("u:", random_term)]]))
    
    # BLUEs
    blue_mat <- matrix(data = fixef, nrow = nrow(blup_mat), ncol = nrow(fixef), byrow = TRUE)
    
    preds <- as.data.frame(blue_mat + blup_mat)
    preds <- cbind(term = row.names(preds), preds)
    names(preds) <- c(random_term, paste0("predicted.", names(preds)[-1]))
    row.names(preds) <- NULL
    
    
  }
  
  # Merge preds with .test; calculate accuracy
  .test1 <- merge(x = .test, y = preds)
  # acc_cor <- cor(.test1[[as.character(fixed)[2]]], .test1$predicted.value)
  # acc_rmse <- sqrt(mean((.test1[[as.character(fixed)[2]]] - .test1$predicted.value)^2))
  # 
  # Create output
  
  out <- list(predictions = .test1)

  return(out)
   
}







# Function that generates regression values for a population against an environmental index
# (i.e. stability), then performs an F test to test differences among the regression coefficients
fwr <- function(data, pheno, gen, env, covariate, weights, model.type = c("indiv", "joint"),
                genotype.random = TRUE, stability.random = FALSE) {
  
  stopifnot(is.data.frame(data))
  cols <- c(pheno, gen, env, covariate, weights)
  col_match <- sapply(cols, `%in%`, names(data))
  if (any(!col_match)) {
    stop("The variables ", paste0(names(col_match)[!col_match], collapse = ", "), " are not columns in 'data'.")
  }
  stopifnot(is.logical(genotype.random))
  stopifnot(is.logical(stability.random))
  
  model.type <- match.arg(model.type)
  
  # Factorize
  data[[gen]] <- as.factor(data[[gen]])
  data[[env]] <- as.factor(data[[env]])
  
  # Drop levels
  data <- droplevels(data)
  
  # Rename weights
  data[["weights"]] <- data[[weights]]
  
  # Unique covariates
  data_cov <- unique(data[c(env, covariate)])
  h_scale <- scale(x = data_cov[[covariate]], center = TRUE, scale = FALSE)

  
  ## Split by a joint model or separate model
  if (model.type == "indiv") {
    
    warning("Note: the 'genotype.random' and 'stability.random' arguments are ignored when model.type = 'indiv'.")
    
    # Create a formula for the individual models
    formula <- reformulate(termlabels = covariate, response = pheno)
    
    # Split by genotype
    data_split <- split(data, data[[gen]])
    
    # Iterate over splits
    output <- list()
    
    for (i in seq_along(data_split)) {
      dat <- data_split[[i]]
      fit <- lm(formula = formula, data = dat)
      output[[i]] <- tibble(term = names(data_split)[i], g = coef(fit)[1], b = coef(fit)[2], d = sum(resid(fit)^2), var_d = sigma(fit)^2)
    }
    
    out <- do.call("rbind", output)
    names(out)[1] <- gen

    # Calculate degrees of freedom
    df_b <- n_distinct(data[[gen]]) - 1
    df_d <- df_b * (n_distinct(data[[env]]) - 2)
    
    # Environmental index sums of squares
    index <- distinct_at(data, vars(all_of(c(env, as.character(formula)[3]))))[[as.character(formula)[3]]]
    e_ss <- sum(index^2)
    
    ## Significance test
    b_ss <- sum(scale(out$b, scale = FALSE)^2) * e_ss
    d_ss <- sum(out$d)
    
    out <- out[c(gen, "g", "b", "var_d")]
    
    # Mean squares
    b_ms <- b_ss / df_b
    d_ms <- d_ss / df_d
    
    # F stat
    Fstat <- b_ms / d_ms
    pvalue <- pf(q = Fstat, df1 = df_b, df2 = df_d, lower.tail = FALSE)
    
    ## Output the results
    tibble(
      regression = list(out),
      Fstat = Fstat,
      pvalue = pvalue
    )
    
  } else {
    
    # require(asreml)
    # require(metafor)
    require(blme)
    require(lmerTest)
    
    # # Fixed and random formula
    # fit_asr <- asreml(fixed = value ~ 1, 
    #                   random = ~ line_name + line_name:environment, data = data, family = asr_gaussian(dispersion = 1.0),
    #                   weights = weights)
    # coef <- summary(fit_asr, coef = TRUE)
    # coef$coef.random[230:250,]
    # fit_asr1 <- asreml(fixed = value ~ 1, 
    #                   random = ~ line_name + line_name:h + line_name:environment, data = data, family = asr_gaussian(dispersion = 1.0),
    #                   weights = weights)
    # coef1 <- summary(fit_asr1, coef = TRUE)
    # coef1$coef.random[460:480,]
    # 
    # # Fit the joint model using metafor
    # mf <- model.frame(reformulate(termlabels = c(gen, covariate), response = pheno), data = data1, weights = weights)
    # y <- model.response(mf)
    # V <- diag(1 / model.weights(mf))
    # mods <- reformulate(termlabels = paste0(gen, ":", covariate))
    # random <- reformulate(termlabels = paste0("1|", gen))
    # 
    # fit_joint_rma <- rma.mv(yi = y, V = V, mods = mods, random = random, data = mf)
    
    
    # Create the formula
    genotype_termlabel <- if (genotype.random) paste0("(1|", gen, ")") else gen
    stability_termlabel <- if (stability.random) paste0("(0 + ", covariate, "|", gen, ")") else paste0(gen, ":", covariate)
    resid_termlabel <- paste0("(1|", gen, ":", env, ")")
    
    formula <- reformulate(termlabels = c(genotype_termlabel, stability_termlabel, resid_termlabel), response = pheno)
    
    # Fit the model
    suppressMessages(fit_blmer <- blmer(formula = formula, data = data, weights = weights, resid.prior = point(1.0), cov.prior = NULL))
    varcor_df <- as.data.frame(VarCorr(fit_blmer))
    # Scale the b variance by the sum of squares of h
    varcor_df$vcov1 <- varcor_df$vcov
    varcor_df$vcov1[which(varcor_df$var1 == "h")] <- varcor_df$vcov[which(varcor_df$var1 == "h")] * mean(h_scale^2)
    # Add the average stage 1 residuals variance
    varcor_df$vcov1[which(varcor_df$grp == "Residual")] <- mean(1 / data[[weights]])
    # # Calculate PVE
    # varcor_df$pve <- varcor_df$vcov1 / sum(varcor_df$vcov1)
    
    
    # Run ranova
    suppressWarnings(fit_blmer_ranova <- ranova(model = fit_blmer))
    fit_blmer_ranova_df <- cbind(term = row.names(fit_blmer_ranova), as_tibble(fit_blmer_ranova))

    # Get the random effects
    randomeff <- ranef(fit_blmer)
    # Get the fixed effects
    fixedeff <- fixef(fit_blmer)
    
    # Get the genotype means
    if (genotype.random) {
      g <- randomeff[[gen]]$`(Intercept)`
      # Get PEV
      g_pev <- attr(randomeff[[gen]], "postVar")$`(Intercept)`[,,]
      # Get variance component
      varG <- subset(varcor_df, grp == paste0(gen, ".1"), vcov, drop = TRUE)
      # Calculate reliability
      g_r2 <- 1 - (g_pev / varG)
      
      
    } else {
      g <- fixedeff[grepl(pattern = gen, x = names(fixedeff))]
      g_r2 <- NULL
      
    }
    
    # Get the stability estimates
    if (stability.random) {
      b <- randomeff[[gen]][[covariate]]
      # Get PEV
      b_pev <- attr(randomeff[[gen]], "postVar")[[covariate]][,,]
      # Get variance component
      varB <- subset(varcor_df, grp == gen, vcov, drop = TRUE)
      # Calculate reliability
      b_r2 <- 1 - (b_pev / varB)
      
    } else {
      b <- fixedeff[grepl(pattern = paste0(":", covariate), x = names(fixedeff))]
      b_r2 <- NULL
      
    }
    
    # Calculate deviation stability
    delta <- randomeff[[paste0(gen, ":", env)]]$`(Intercept)`
    delta_pev <- attr(randomeff[[paste0(gen, ":", env)]], "postVar")[,,]
    # Get variance component
    varDelta <- subset(varcor_df, grp == paste0(gen, ".", env), vcov, drop = TRUE)
    # Calculate reliability
    delta_r2 <- 1 - (delta_pev / varDelta)
    
    data1 <- data[gen]
    data1$delta <- delta
    data1$delta_dereg <- delta / delta_r2
    data1$delta_r2 <- delta_r2
    
    if (all(data1$delta == 0)) {
      out <- distinct(data1)
      names(out)[-1] <- c("var_d", "var_d_dereg", "delta_r2")
    } else {
      out <- aggregate(reformulate(termlabels = gen, response = "cbind(delta, delta_dereg)"), data = data1, FUN = var)
      names(out)[-1] <- c("var_d", "var_d_dereg")
      out2 <- aggregate(reformulate(termlabels = gen, response = "delta_r2"), data = data1, FUN = mean)
      out <- merge(out, out2)
    }
    
    out$g <- g + fixedeff[1]
    out$g_r2 <- g_r2
    out$b <- b
    out$b_r2 <- b_r2
    
    out <- as_tibble(out)
    
    tibble(
      regression = list(out),
      vars = list(varcor_df),
      ranova = list(fit_blmer_ranova_df)
    )

  }
  
  
}





prune_LD2 <- function(x, r2.max = 0.80) {
  
  ## Error checking
  # Check the marker matrix
  # stopifnot(check_marker_matrix(x))
  # check r2max
  stopifnot(r2.max >= 0 & r2.max <= 1)
  
  # calculate minor allele frequency
  maf <- calc_maf(x)
  
  # calculate the correlation between all markers; square it
  if (any(is.na(x))) {
    all_marker_r <- cor(x, use = "pairwise.complete.obs")^2
    
  } else {
    all_marker_r <- cor(x)^2
  }
  
  # Set the lower half (including diagonal) to NA
  all_marker_r1 <- all_marker_r
  all_marker_r1[lower.tri(all_marker_r1, diag = TRUE)] <- NA
  
  # Get a matrix of those entries that are elevated
  elevated_marker_r <- all_marker_r1 > r2.max
  # Get the coordinates of those entries
  which_elevated <- which(x = elevated_marker_r, arr.ind = TRUE)
  
  markers_remove <- character()
  # While loop
  i = 1
  while(nrow(which_elevated) > 0) {
    
    # Subset the first row
    coords <- which_elevated[1,]
    # Extract the coordinate
    r2_coord <- all_marker_r1[coords[1], coords[2], drop = FALSE]
    
    # marker pair
    markers <- unlist(dimnames(r2_coord))
    # Identify the marker with higher MAF
    higher_maf <- which.max(maf[markers])
    
    # Toss that marker
    marker_remove <- names(higher_maf)
    markers_remove[i] <- marker_remove
    
    # Find the row/col containing this marker
    row_remove <- col_remove <- which(row.names(all_marker_r1) == marker_remove)
    
    which_elevated <- subset.matrix(which_elevated, which_elevated[,"row"] != row_remove & which_elevated[,"col"] != col_remove)
    
    # advance i
    i <- i + 1
    
  }
  
  # Remove the markers from the marker matrix
  cols_keep <- setdiff(seq_len(ncol(x)), which(colnames(x) %in% markers_remove))
  x[,cols_keep,drop = FALSE]
  
}


calc_maf2 <- function(x) {
  
  ## Error checking
  # Check the marker matrix
  # stopifnot(check_marker_matrix(x))
  
  # Filter for MAF
  af <- colMeans(x + 1, na.rm = TRUE) / 2
  maf <- pmin(af, 1 - af)
  
  return(maf)
  
}


filter_snps2 <- function(x, r2.max, maf.min, indiv.miss.max, snp.miss.max) {
  
  ## Error checking
  # Check the marker matrix
  # stopifnot(check_marker_matrix(x))
  # check r2max
  
  ## Filter on missingness
  if (!missing(snp.miss.max)) {
    stopifnot(snp.miss.max >= 0 & snp.miss.max <= 1)
    
    snp_missing <- colMeans(is.na(x))
    x1 <- x[,snp_missing <= snp.miss.max, drop = FALSE]
    
  } else {
    x1 <- x
  }
  
  if (!missing(indiv.miss.max)) {
    stopifnot(indiv.miss.max >= 0 & indiv.miss.max <= 1)
    
    indiv_missing <- rowMeans(is.na(x1))
    x2 <- x1[indiv_missing <= indiv.miss.max, , drop = FALSE]
    
  } else {
    x2 <- x1
  }
  
  ## Filter on maf
  if (!missing(maf.min)) {
    stopifnot(maf.min >= 0 & maf.min <= 1)
    
    # calculate minor allele frequency
    maf <- calc_maf2(x2)
    x3 <- x2[,maf >= maf.min, drop = FALSE]
    
  } else {
    x3 <- x2
  }
  
  if (!missing(r2.max)) {
    stopifnot(r2.max >= 0 & r2.max <= 1)
    # Prune on LD
    x4 <- prune_LD2(x = x3, r2.max = r2.max)
    
  } else {
    x4 <- x3
  }
  
  # Return x4
  return(x4)
  
}

















