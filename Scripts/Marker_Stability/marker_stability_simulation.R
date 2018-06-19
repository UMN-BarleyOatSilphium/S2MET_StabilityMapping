## S2MET Mapping
## 
## Simulation of marker effect stability
## This script will run a simulation to explore the additive model of marker effect
## reaction norms towards phenotypic reaction norms
## 
## 
## Author: Jeff Neyhart
## Last updated: June 17, 2018
## 

# Run the source script - local
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

library(qtl)
library(pbsim)

# Get the barley data
load("C:/Users/Jeff/GoogleDrive/BarleyLab/Projects/Side Projects/Resources/s2_cap_simulation_data.RData")

n_iter <- 50
pop_size <- 250
n_env <- 5
n_rep <- 1
nqtl <- 100
# Trait heritability
h2 <- 0.5

## Use the cap population
### Create the genome and genetic model
gen_map <- s2_snp_info %>% 
  select(1,3,4) %>% 
  data.frame(row.names = .$rs) %>% 
  select(-1) %>% 
  table_to_map() %>% 
  jittermap()

# Select individuals to retain
pop_ind <- s2_cap_genos[sort(sample(row.names(s2_cap_genos), pop_size)),]

## Sample 100 QTL from the 6th and 7th chromosomes
gen_map1 <- c(gen_map[1:5], map(gen_map[6:7], ~sort(sample(., nqtl))))

# Create the genome
genome <- sim_genome(map = gen_map1)

## Replicate
sim_out <- replicate(n = n_iter, expr = {

  # Simulate a genetic model
  qtl_model <- matrix(NA, ncol = 4, nrow = nqtl)
  qtl_model1 <- cbind(chr = 6, pos = gen_map1$`6H`, add_eff = rnorm(nqtl), dom_eff = 0)
  qtl_model2 <- cbind(chr = 7, pos = gen_map1$`7H`, add_eff = rnorm(nqtl), dom_eff = 0)
  
  genome1_1 <- sim_gen_model(genome = genome, qtl.model = qtl_model, add.dist = "normal")
  genome1_2 <- sim_gen_model(genome = genome, qtl.model = qtl_model, add.dist = "normal")
  
  genome1 <- genome1_1
  genome1$gen_model <- list(genome1_1$gen_model[[1]], genome1_2$gen_model[[1]])
  genome1$gen_model[[1]]$pos <- gen_map1$`6H`
  genome1$gen_model[[1]]$qtl_name <- names(gen_map1$`6H`)
  genome1$gen_model[[1]]$chr <- 6
  
  genome1$gen_model[[2]]$pos <- gen_map1$`7H`
  genome1$gen_model[[2]]$qtl_name <- names(gen_map1$`7H`)
  genome1$gen_model[[2]]$chr <- 7

  # Create the population
  pop <- create_pop(genome = genome1, geno = pop_ind[,markernames(genome1, include.qtl = TRUE)])
  # Get the genotypic values and standardize - these are for stability only
  g <- pop$geno_val$trait1
  
  varG <- var(g)
  # Intended residual variance
  varR_tot <- (varG / h2) - varG
  varGE <- varR_tot * 0.1
  varR <- varR_tot - varGE
  # Intended environmental variance
  varE <- varG * 8
  
  # Remember dividing values by a number reduces the variance by the square
  # of that number
  b <- c(scale(pop$geno_val$trait2, scale = sqrt(var(b) / varGE)))
  
  # Genotype
  M <- genotype(genome = genome1, pop = pop)
  

  
  # # Generate reaction norm values
  # b <- rnorm(n = pop_size, mean = 0, sd = sqrt(varGE))
  # Generate environmental effects
  h <- sort(rnorm(n = n_env, mean = 0, sd = sqrt(varE)))
  # # Generate genotypic effects (slopes)
  # g <- rnorm(n = pop_size, mean = 0, sd = sqrt(varG))
  
  # Generate fitted genotypic values in each environment
  bh <- outer(b, h)
  # Generate random deviations from the regression line
  ep <- replicate(n = n_env * n_rep, rnorm(n = pop_size, mean = 0, sd = sqrt(varR)))
  
  y <- 0 + replicate(n = n_env * n_rep, g) + 
    matrix(h, ncol = n_rep * n_env, nrow = pop_size, byrow = T) + 
    matrix(bh, ncol = n_rep * n_env, nrow = pop_size, byrow = F) +
    ep
    
  # Convert to df
  dimnames(y) <- list(
    indnames(pop),
    paste(rep(paste0("env", seq(n_env)), n_rep), rep(paste0("rep", seq(n_rep)), each = n_env), sep = "_")
  )
  
  y_df <- as.data.frame(y) %>% 
    rownames_to_column("gen") %>% 
    gather(obs, y, -gen) %>% 
    separate(obs, c("env", "rep"), sep = "_")
  
  fit <- lm(y ~ env + gen, y_df, contrasts = list(gen = "contr.sum", env = "contr.sum"))
  beta <- coefficients(fit)[2:n_env]
  h_hat <- c(beta, -sum(beta))
  
  # Calculate reaction norms
  y_df1 <- y_df %>% 
    group_by(gen, env) %>%
    summarize(y = mean(y))
  
  phen_reac_norms <- y_df1 %>% 
    do(b_hat = coef(lm(y ~ h, data = .))[2]) %>% 
    unnest()
  
  # Calculate marker specific reaction norms
  marker_env_eff <- y_df1 %>% 
    group_by(env) %>% 
    do(data.frame(marker = colnames(M), u = mixed.solve(y = .$y, Z = M)$u)) %>%
    mutate(marker = factor(marker, levels = colnames(M)))
  
  marker_reac_norms <- marker_env_eff %>% 
    group_by(marker) %>% 
    do(b_hat = coef(lm(u ~ h, data = .))[2]) %>% 
    unnest()
  
  ## Now use the marker reaction norms as effect of phenotypic reaction norms
  b_pred <- M %*% marker_reac_norms$b_hat
  # Correlate
  b_cor <- cor(phen_reac_norms$b_hat, b_pred)
  
  ## Now calculate b marker effects
  b_marker_effs <- mixed.solve(y = phen_reac_norms$b_hat, Z = M)$u
  # Correlate
  b_marker_cor <- cor(marker_reac_norms$b_hat, b_marker_effs)
  
  
  c(b_cor = b_cor, b_marker_cor = b_marker_cor)
  
})


#
