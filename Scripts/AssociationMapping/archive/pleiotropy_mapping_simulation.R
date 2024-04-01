## Pleiotropy GWAS simulation

# Run the source script - local
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

library(pbsim)

n_iter <- 50

# The below simulations estimates the empirical FDR when the NULL hypothesis of 
# no QTL is TRUE. We do this by putting all the QTL on one chromosome, and testing
# only markers on another.
empirical_fdr <- replicate(n = n_iter, expr = {

  ## Generate a genome of 1 chromosome and 100 markers
  ## All of the QTL will be on chromosome 1, and all the markers on chromosome 2
  genome <- sim_genome(len = rep(150, 1), n.mar = c(100))
  # Simulate QTL positions
  genome1 <- sim_multi_gen_model(genome = genome, add.dist = "normal", corr = 0.5, prob.corr = cbind(0, 1),
                                 qtl.model = replicate(2, matrix(NA, ncol = 4, nrow = 100), simplify = FALSE))
  
  # Add a second chromosome
  genome1$len <- c(genome1$len, genome1$len)
  genome1$n.mar <- c(genome1$n.mar, 2000)
  genome1$map <- list(`1` = genome1$map$`1`, `2` = qtl::sim.map(len = c(150, 150), n.mar = 2000, include.x = FALSE)[[2]])
  
  # Simulate the population
  pop <- sim_pop(genome = genome1, n.ind = 200)
  
  # Phenotype
  pop1 <- sim_phenoval(pop = pop, h2 = 0.7, n.env = 2, n.rep = 2)
  
  # Genotype
  pop_geno <- genotype(genome = genome1, pop = pop1)
  K <- A.mat(X = pop_geno, min.MAF = 0, max.missing = 1)
  
  # Format
  pheno <- pop1$pheno_val$pheno_mean
  geno <- genome1$map$`2` %>% 
    structure(class = "numeric") %>%
    data.frame(marker = names(.), chrom = 2, pos = ., t(pop_geno))
  
  # GWAS
  gwas_mv_out <- mv_gwas(pheno = pheno, geno = geno, K = K, n.PC = 0)
  
  alpha <- 0.05
  
  # Get a matrix of pvalues
  pvalue_mat <- 10^-gwas_mv_out[,-c(1:3)]
  # Get the minimum p value per marker
  max_p <- apply(X = pvalue_mat, MARGIN = 1, FUN = max)
  
  # How many pvalues are both below the significance threshold?
  gwas_mv_reject <- gwas_mv_out %>% 
    mutate_at(vars(contains("trait")), ~10^-. <= alpha) %>%
    filter(trait1 & trait2)
  
  n_false_disc <- nrow(gwas_mv_reject)
  (p_false_disc <- n_false_disc / nrow(gwas_mv_out))
  
})



