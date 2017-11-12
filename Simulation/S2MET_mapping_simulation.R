## S2MET Mapping Simulation
## 
## Run a simulation to test different association mapping models to detect QxE
## 

# Packages and directories
# List of packages
packages <- c("dplyr", "purrr", "tibble", "tidyr", "readr", "stringr", "readxl", "modelr", 
              "parallel", "purrrlyr", "rrBLUP", "pbr", "qtl", "pbsim")

# Set the directory of the R packages
package_dir <- NULL
# package_dir <- "/panfs/roc/groups/6/smithkp/neyha001/R/x86_64-pc-linux-gnu-library/3.4/"

# Load all packages
invisible(lapply(packages, library, character.only = TRUE, lib.loc = package_dir))


## Directories
proj_dir <- "C:/Users/Jeff/Google Drive/Barley Lab/Projects/S2MET_Mapping//"
# proj_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET_Mapping/" 

alt_proj_dir <- "C:/Users/Jeff/Google Drive/Barley Lab/Projects/S2MET"
# alt_proj_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET/"

# Other directories
fig_dir <- file.path(proj_dir, "Figures/")
map_dir <- file.path(proj_dir, "Mapping")
entry_dir <- file.path(alt_proj_dir, "Plant_Materials")
analysis_dir <- file.path(proj_dir, "Analysis")
result_dir <- file.path(proj_dir, "Results")

data("s2_cap_genos")
data("s2_snp_info")
data("s2_cap_line_info")

s2_cap_line_info1 <- s2_cap_line_info %>%
  filter(CAP_Name %in% row.names(s2_cap_genos))

s2_cap_genos1 <- s2_cap_genos[s2_cap_line_info1$CAP_Name,]

## Summarize the line info
program_info <- s2_cap_line_info1 %>% 
  select(line = CAP_Name) %>% 
  mutate(program = str_extract(line, "[A-Z]{1,}[0-9]{0,1}")) %>%
  group_by(program) %>%
  mutate(proportion = n() / nrow(.)) %>%
  ungroup()


### Create the genome and genetic model
gen_map <- s2_snp_info %>% 
  select(1,3,4) %>% 
  data.frame(row.names = .$rs) %>% 
  select(-1) %>% 
  table_to_map() %>% 
  jittermap()

genome <- sim_genome(map = gen_map)

# Pull out the map as a df
use_map <- map_to_table(genome) %>%
  rownames_to_column("marker")



### Set some simulation parameters
n_qtl_list <- c(30, 100) # Number of QTL
h2_list <- c(0.5, 0.8) # Heritability
n_pop_list <- c(175, 350, 700, 1400)

n_iter <- 10 # Number of simulation iterations


### Contstant parameters
n_env <- 2 # Number of environments
n_rep <- 1 # Number of reps
p_qxe <- 0.50 # Proportion QTL with QxE
varGE_scale <- 0.5 # Scaling parameter of varGE to varR
models <- c("K", "G", "KE", "GE") # Association models
sig_cutoff <- c(0.01, 0.05, 0.10) # Significance thresholds

## Data frame of parameters
param_df <- expand.grid(n_qtl = n_qtl_list, h2 = h2_list, n_pop = n_pop_list)


# Detect cores
n_cores <- detectCores()

# Iterate over the parameter df and run on cores
sim_results <- param_df %>%
  mutate(core = sort(rep(seq(n_cores), length.out = nrow(.)))) %>%
  split(.$core) %>%
  # Apply over cores
  mclapply(X = ., FUN = function(trial) {

    n_qtl <- trial$n_qtl
    h2 <- trial$h2
    n_pop <- trial$n_pop
    
    # Create list
    trial_list <- vector("list", n_iter)
    
    # Iterate over replications
    for (iter in seq(n_iter)) {
    
      ## Randomly sample the S2 population weighted by the proportion of each sub-population
      ## in the large population
      sample_pop <- program_info %>% 
        sample_n(tbl = ., size = n_pop, weight = proportion) %>%
        arrange(line)
      
      pop_geno <- s2_cap_genos1[sample_pop$line,]
        
        
      ### Simulate a genetic model
      genome1 <- sim_gen_model(genome = genome, qtl.model = matrix(NA, nrow = n_qtl, ncol = 4),
                               add.dist = "geometric")
      
      # Create a population
      assoc_panel <- create_pop(genome = genome1, geno = pop_geno)
      
      # Randomly select p_qxe of the QTL
      sample_qxe <- pull_qtl(genome1) %>% 
        sample_frac(size = p_qxe)
      
      
      ### Phenotype
      
      # Get the genetic variance
      varG <- var(assoc_panel$geno_val$trait1)
      # Scale the residual variance
      varR <- (((varG / h2) - varG) * n_env) / (1 + varGE_scale)
      varGE <- varR * varGE_scale
      
      # Should varGE be divided by the number of effective QTL?
      
      # Genotypic values
      geno_value <- matrix(data = assoc_panel$geno_val$trait1, nrow = n_pop, ncol = n_env * n_rep)
      
      # Assigne environment-specific QTL effects for the QTLxE
      sample_qxe_effects <- sample_qxe %>% 
        by_row(.d = ., ..f = function(i) rnorm(n = n_env, mean = 0, sd = sqrt(varGE / (p_qxe * n_qtl))), 
               .collate = "cols", .to = "env") %>% 
        left_join(pull_qtl(genome1), ., by = c("chr", "pos", "add_eff", "dom_eff", "qtl_name", "qtl1_pair", "trait")) %>%
        mutate_at(vars(starts_with("env")), ~ ifelse(is.na(.), 0, .))
      
      # Pull out QTL genotypes
      qtl_geno <- pull_genotype(genome = genome1, geno = pop_geno, loci = sample_qxe_effects$qtl_name) - 1
      
      # Calculate the environment-specific genotype values
      gxe_value <- qtl_geno %*% as.matrix(select(sample_qxe_effects, starts_with("env")))
      gxe_value <- matrix(data = gxe_value, ncol = n_env * n_rep, nrow = nrow(gxe_value))
      
      # Simulate environmental effects
      e_value <- matrix(data = rnorm(n = n_env, mean = 0, sd = sqrt(varG * 8)), nrow = n_pop, ncol = n_env * n_rep, byrow = T)
      
      # Simulate residual effects
      resid_value <- matrix(data = rnorm(n = n_pop * n_env * n_rep, mean = 0, sd = sqrt(varR)), nrow = n_pop, ncol = n_env * n_rep)
      
      # Add for phenotypes
      pheno_value <- geno_value + e_value + gxe_value + resid_value
      
      ### Modeling
      # phenotypes
      pheno_tomodel <- pheno_value %>% 
        as.data.frame() %>%
        set_names(paste(paste("env", seq(n_env), sep = ""), rep(paste("rep", seq(n_rep), sep = ""), each = n_env), sep = "_")) %>%
        mutate(line = indnames(assoc_panel)) %>%
        gather(obs, pheno, -line) %>%
        separate(obs, c("env", "rep"), sep = "_") %>%
        select(-rep)
      
      
      # ## Run linear model to check variance components
      # fit <- lm(pheno ~ line*env, data = pheno_tomodel)
      # tidy_aov <- tidy(anova(fit))
      # varR_hat <- subset(tidy_aov, term == "Residuals", meansq, drop = T)
      # varGE_hat <- (subset(tidy_aov, term == "line:env", meansq, drop = T) - varR_hat) / n_rep
      # varG_hat <- (subset(tidy_aov, term == "line", meansq, drop = T) - subset(tidy_aov, term == "line:env", meansq, drop = T)) / (n_rep * n_env)
      # 
      
      
      # Genotype the panel
      assoc_panel_geno <- genotype(genome = genome1, pop = assoc_panel) %>%
        t() %>%
        data.frame(rs = row.names(.), ., row.names = NULL, stringsAsFactors = FALSE,
                   check.names = FALSE) %>%
        rename(marker = rs)
      
      ## Combine genotypes with snp information
      geno_tomodel <- right_join(x = use_map, assoc_panel_geno, "marker") %>%
        as.data.frame()
      
      ### Association analysis
      gwas_main <- map(models, ~gwas(pheno = pheno_tomodel, geno = geno_tomodel, 
                                     fixed = ~ env, model = ., n.PC = 2, P3D = TRUE, 
                                     n.core = 1, test.qxe = T)) %>% set_names(models)
      
      # Adjust p_values for main effect
      gwas_main_adj_p <- gwas_main %>%
        map("scores") %>% 
        map(~filter(., term == "main_effect") %>% 
              mutate(adj_p_value = p.adjust(p_value, "fdr")))
      # Now for QxE
      gwas_qxe_adj_p <- gwas_main %>%
        map("scores") %>% 
        map(~filter(., term == "qxe") %>% 
              mutate(adj_p_value = p.adjust(p_value, "fdr")))
      
      # Significant markers - map across different significance levels
      sig_mar_main <- map(sig_cutoff, function(sig) 
        map(gwas_main_adj_p, ~filter(., adj_p_value <= sig))) %>% 
        set_names(paste("fdr", sig_cutoff, sep = ""))
      
      sig_mar_qxe <- map(sig_cutoff, function(sig) 
        map(gwas_qxe_adj_p, ~filter(., adj_p_value <= sig))) %>% 
        set_names(paste("fdr", sig_cutoff, sep = ""))
      
      
      # Pull out the map and mark loci as a marker or a QTL
      loci <- use_map %>%
        left_join(., select(sample_qxe_effects, marker = qtl_name, chr = chr, pos, add_eff, starts_with("env")),
                  by = c("marker", "chr", "pos"))
      
      # Add significance to the mapped markers
      loci_main <- sig_mar_main %>%
        map(., ~map(., select, marker, chr = chrom, pos, p_value) %>% 
              map(~left_join(loci, ., by = c("marker", "chr", "pos"))))
      
      loci_qxe <- sig_mar_qxe %>%
        map(., ~map(., select, marker, chr = chrom, pos, p_value) %>% 
              map(~left_join(loci, ., by = c("marker", "chr", "pos"))))
      
      # Assign as a marker or QTL and significant
      loci1_main <- loci_main %>% 
        map(., ~map(., ~mutate(., type = if_else(is.na(add_eff), "marker", "qtl"),
                    sig = !is.na(p_value))))
      
      loci1_qxe <- loci_qxe %>% 
        map(., ~map(., ~mutate(., type = if_else(is.na(add_eff), "marker", "qtl"),
                               sig = !is.na(p_value))))
      
      
      ## Find the false positives for the main effect
      # If a marker is significant, but there is no adjacent QTL, declare a false positive
      # If a marker is not significant, but there is no adjacent QTL, declare a true negative
      # If a marker is signficant, but there is an adjacent QTL, declare a true positive
      # If a marker is not significant, but there is an adjacent marker, declare a false negative
      
      # Split by chromosome
      mar_assoc_main <- loci1_main %>%
        map(., ~map(., function(loc) {
          loc %>%
            split(.$chr) %>% 
            map(function(loci_map) {
              
              assoc_results <- character(nrow(loci_map))
              for (i in seq_along(assoc_results)) {
                if (loci_map[i,"type"] == "qtl") {
                  assoc_results[i] <- NA
                  
                } else if (loci_map[i,"sig"]) {
                  assoc_results[i] <- ifelse(any(loci_map[c(i - 1, i + 1),"type"] == "qtl", na.rm = T), 
                                             "true_positive", "false_positive")
                  
                } else {
                  # Is there a QTL in the adjacent interval?
                  assoc_results[i] <- ifelse(any(loci_map[c(i - 1, i + 1),"type"] == "qtl", na.rm = T), 
                                             "false_negative", "true_negative")
                  
                  }
              }
              data.frame(loci_map, assoc_results = assoc_results)
              }) %>%
            do.call("rbind", .) }))
      
      ## Find the false positives for the qxe effect
      # If a marker is significant, but there is no adjacent QTL, declare a false positive
      # If a marker is not significant, but there is no adjacent QTL, declare a true negative
      # If a marker is signficant, but there is an adjacent QTL, declare a true positive
      # If a marker is not significant, but there is an adjacent marker, declare a false negative
      
      # Split by chromosome
      mar_assoc_qxe <- loci1_qxe %>%
        map(., ~map(., function(loc) {
          loc %>%
            split(.$chr) %>% 
            map(function(loci_map) {
              
              assoc_results <- character(nrow(loci_map))
              for (i in seq_along(assoc_results)) {
                if (loci_map[i,"type"] == "qtl") {
                  assoc_results[i] <- NA
                  
                  # A marker is a true positive if it is significant, there is a QTL
                  # in the adjacent interval, and the QTL has a non-zero QxE effect
                } else if (loci_map[i,"sig"]) {
                  # Is there a QTL?
                  is_qtl <- loci_map[c(i - 1, i + 1),"type"] == "qtl"
                  # Is there qxe effect
                  is_qxe <- rowSums(loci_map[c(i - 1, i + 1), c("env1", "env2")] != 0) > 0
                  
                  assoc_results[i] <- ifelse(any(is_qtl & is_qxe, na.rm = T), "true_positive", "false_positive")
                  
                } else {
                  # Is there a QTL?
                  is_qtl <- loci_map[c(i - 1, i + 1),"type"] == "qtl"
                  # Is there qxe effect
                  is_qxe <- rowSums(loci_map[c(i - 1, i + 1), c("env1", "env2")] != 0) > 0
                  
                  # Is there a QTL in the adjacent interval?
                  assoc_results[i] <- ifelse(any(is_qtl & is_qxe, na.rm = T), "false_negative", "true_negative")
                  
                }
              }
              data.frame(loci_map, assoc_results = assoc_results)
            }) %>%
            do.call("rbind", .) }))
      
      
      
      ## Determine the number of true QTL detected
      # If a QTL is flanked by at least one significant marker, the QTL is discovered
      # Split by chromosome
      qtl_discover_main <- loci1_main %>%
        map(., ~map(., function(loc) {
          loc %>%
            split(.$chr) %>% 
            map(function(loci_map) {
              
              # Pull out markers and marker interval
              mar_map <- subset(x = loci_map, type == "marker")
              mar_interval <- mar_map$pos
              
              # Determine if QTL were discovered
              subset(x = loci_map, type == "qtl") %>% 
                mutate(first_mar = findInterval(x = pos, vec = mar_interval), 
                       sec_mar = first_mar + 1, discovered = mar_map[first_mar,"sig"] | mar_map[sec_mar,"sig"]) %>%
                select(marker:env2, discovered) }) %>%
            bind_rows() }))
      
      qtl_discover_qxe <- loci1_qxe %>%
        map(., ~map(., function(loc) {
          loc %>%
            split(.$chr) %>% 
            map(function(loci_map) {
              
              # Pull out markers and marker interval
              mar_map <- subset(x = loci_map, type == "marker")
              mar_interval <- mar_map$pos
              
              # Determine if QTL were discovered
              filter(loci_map, type == "qtl", env1 != 0 | env2 != 0) %>% 
                mutate(first_mar = findInterval(x = pos, vec = mar_interval), 
                       sec_mar = first_mar + 1, discovered = mar_map[first_mar,"sig"] | mar_map[sec_mar,"sig"]) %>%
                select(marker:env2, discovered) }) %>%
            bind_rows() }))
      
      # Separate true/false and negative/positive
      mar_assoc_main1 <- mar_assoc_main %>% 
        as_data_frame() %>%
        mutate(model = models, term = "main_effect") %>% 
        gather(sig_cutoff, data, -model, -term)
      
      mar_assoc_qxe1 <- mar_assoc_qxe %>% 
        as_data_frame() %>%
        mutate(model = models, term = "qxe") %>% 
        gather(sig_cutoff, data, -model, -term)
      
      assoc_count <- bind_rows(mar_assoc_main1, mar_assoc_qxe1) %>%
        unnest() %>%
        separate(col = assoc_results, into = c("true_false", "neg_pos"), sep = "_") %>% 
        filter(type == "marker") %>% 
        group_by(model, sig_cutoff, term, true_false, neg_pos) %>% 
        summarize(n = n())
      
      # Count the number of QTL that were discovered
      discovered_count <- map(list(qtl_discover_main, qtl_discover_qxe), ~as_data_frame(.) %>%
                                mutate(model = models) %>% gather(sig_cutoff, data, -model) %>% 
                                unnest() %>% group_by(model, sig_cutoff) %>% 
                                summarize(n_discovered = sum(discovered), prop_discovered = mean(discovered))) %>%
        set_names(c("main_effect", "qxe"))
      
      # Add DF to the list
      trial_list[[iter]] <- data_frame(iter = iter, marker_association = list(assoc_count), qtl_discovery = list(discovered_count))
      
    } # Iteration loop close
    
    # Data.frame of model and list of contingency data.frames
    trial_list %>% 
      bind_rows() %>% 
      gather(datatype, data, -iter)
    
  }, mc.cores = n_cores)
    
      
# Save the results
save_file <- file.path(result_dir, "S2MET_mapping_gwas_qxe_simulation.RData")


# ## Manipulate data
# test <-trial_list_df %>% 
#   filter(datatype == "marker_association") %>% 
#   unnest() %>% select(-datatype) %>% 
#   complete(iter, model, sig_cutoff, term, true_false, neg_pos) %>% 
#   filter(neg_pos == "positive") %>% 
#   spread(true_false, n) %>% 
#   group_by(model, sig_cutoff, term) %>% 
#   summarize(n_pos_mar = mean(true, na.rm = T), n_pos_mar_sd = sd(true, na.rm = T),
#             n_false_pos = mean(false, na.rm = T), n_false_pos_sd = sd(false, na.rm = T))
# 
