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
package_dir <- "/panfs/roc/groups/6/smithkp/neyha001/R/x86_64-pc-linux-gnu-library/3.4/"

# Load all packages
invisible(lapply(packages, library, character.only = TRUE, lib.loc = package_dir))


## Directories
proj_dir <- "C:/Users/Jeff/Google Drive/Barley Lab/Projects/S2MET_Mapping//"
proj_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET_Mapping/" 

alt_proj_dir <- "C:/Users/Jeff/Google Drive/Barley Lab/Projects/S2MET"
alt_proj_dir <- "/panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET/"

# Other directories
fig_dir <- file.path(proj_dir, "Figures/")
map_dir <- file.path(proj_dir, "Mapping")
entry_dir <- file.path(alt_proj_dir, "Plant_Materials")
analysis_dir <- file.path(proj_dir, "Analysis")
result_dir <- file.path(proj_dir, "Results")
sim_dir <- file.path(proj_dir, "Simulation")

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
n_pop_list <- c(175, 350, 700)

n_iter <- 100 # Number of simulation iterations
max_qtl <- max(n_qtl_list)


### Contstant parameters
n_env <- 2 # Number of environments
n_rep <- 1 # Number of reps
p_qxe <- 0.50 # Proportion QTL with QxE
varGE_scale <- 0.5 # Scaling parameter of varGE to varR
models <- c("K", "G", "KE", "GE", "KKE", "GGE") # Association models
sig_cutoff <- c(0.01, 0.05, 0.10) # Significance thresholds

## Data frame of parameters
param_df <- expand.grid(iter = seq(n_iter), n_qtl = n_qtl_list, h2 = h2_list, n_pop = n_pop_list)


## Subset 
param_df <- param_df %>%
  filter( (n_qtl == 30 & h2 == 0.8) | (n_qtl == 100 & h2 == 0.5) )

# Detect cores
n_cores <- detectCores()

# Iterate over the parameter df and run on cores
sim_results <- param_df %>%
  mutate(core = sort(rep(seq(n_cores), length.out = nrow(.)))) %>%
  split(.$core) %>%
  # Apply over cores
  mclapply(X = ., FUN = function(trial) {
    # Apply down rows
    by_row(trial, function(row) {

      n_qtl <- row$n_qtl
      h2 <- row$h2
      n_pop <- row$n_pop
      iter <- row$iter
      
      
      ## Randomly sample the S2 population weighted by the proportion of each sub-population
      ## in the large population
      sample_pop <- program_info %>% 
        sample_n(tbl = ., size = n_pop, weight = proportion) %>%
        arrange(line)
      
      pop_geno <- s2_cap_genos1[sample_pop$line,]
        
        
      ### Simulate a genetic model
      genome1 <- sim_gen_model(genome = genome, qtl.model = matrix(NA, nrow = n_qtl, ncol = 4),
                               add.dist = "geometric", max.qtl = max_qtl)
      
      # Create a population
      assoc_panel <- create_pop(genome = genome1, geno = pop_geno)
      
      # Randomly select p_qxe of the QTL
      sample_qxe <- pull_qtl(genome1) %>% 
        filter(add_eff != 0) %>%
        sample_frac(size = p_qxe) %>%
        arrange(chr, pos)
      
      
      ### Phenotype
      
      # Get the genetic variance
      varG <- var(assoc_panel$geno_val$trait1)
      # Scale the residual variance
      varR <- (((varG / h2) - varG) * n_env) / (1 + varGE_scale)
      varGE <- varR * varGE_scale
      
    
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
      
      ## Adjust p-values
      gwas_adj_p <- gwas_main %>% 
        map(., ~mutate(.$scores, model = .$metadata$model)) %>% 
        bind_rows() %>%
        group_by(term, model) %>% 
        mutate(adj_p_value = p.adjust(p_value, method = "fdr")) %>%
        ungroup()
      
      # Determine significant markers
      gwas_sig_mar <- gwas_adj_p %>% 
        mutate(sig = list(as.data.frame(t(set_names(sig_cutoff, sig_cutoff))))) %>% 
        unnest() %>% 
        rename_at(vars(starts_with("0")), ~str_c("fdr", .)) %>% 
        mutate_at(vars(starts_with("fdr")), funs(sig = adj_p_value <= .))
      
      ## Split markers by model and by term
      gwas_sig_mar_split <- gwas_sig_mar %>%
        split(list(.$term, .$model))
      
  
      # Pull out the map and mark loci as a marker or a QTL
      loci <- use_map %>%
        left_join(., select(sample_qxe_effects, marker = qtl_name, chr = chr, pos, add_eff, starts_with("env")),
                  by = c("marker", "chr", "pos"))
      
      # Add significance to the mapped markers
      loci_sig <- gwas_sig_mar_split %>% 
        map(., ~left_join(loci, ., by = c("marker", "chr" = "chrom", "pos"))) %>% 
        map(~mutate(., type = if_else(is.na(add_eff), "marker", "qtl"))) %>%
        # Remove the QTL with no effect (i.e. masked)
        map(~filter(., add_eff != 0 & type == "qtl" | type == "marker"))
      
      ## Find the false positives for the main effect
      # If a marker is significant, but there is no adjacent QTL, declare a false positive
      # If a marker is not significant, but there is no adjacent QTL, declare a true negative
      # If a marker is signficant, but there is an adjacent QTL, declare a true positive
      # If a marker is not significant, but there is an adjacent marker, declare a false negative
      
      # Split the loci by main effect or QxE
      loci_sig_main <- loci_sig[str_detect(names(loci_sig), "main_effect")]
      loci_sig_qxe <- loci_sig[str_detect(names(loci_sig), "qxe")]
      
      
      mar_assoc_main <- loci_sig_main %>% 
        map(~split(., .$chr) %>% 
            map(function(loci_map) {
              loci_map <- as.data.frame(loci_map) %>%
                arrange(pos)
              
              assoc_results <- vector("list", nrow(loci_map))
              for (i in seq_along(assoc_results)) {
                if (loci_map[i,"type"] == "qtl") {
                  assoc_results[[i]] <- NA
                  
                } else if (is.na(loci_map[i,"p_value"])) {
                  
                  assoc_results[[i]] <- NA
                  
                } else {
                  
                  # IS there a QTL in the interval?
                  is_qtl <- any(loci_map[c(i - 1, i + 1),"type"] == "qtl", na.rm = T)
                  # What are the sinificance levels
                  sig_levels <- loci_map[i,c("fdr0.01_sig", "fdr0.05_sig", "fdr0.1_sig")]
                  
                  if (is_qtl) {
                    assoc_results[[i]] <- ifelse(sig_levels, "true_positive", "false_negative")
                    
                  } else {
                    assoc_results[[i]] <- ifelse(sig_levels, "false_positive", "true_negative")
                    
                  }
                }
              }
              cbind(select(loci_map, -contains("fdr")), 
                    `dimnames<-`(x = do.call("rbind", assoc_results), list(NULL, str_c("fdr", sig_cutoff))))
            }) %>%
            do.call("rbind", .) )
      
      ## Find the false positives for the qxe effect
      # If a marker is significant, but there is no adjacent QTL, declare a false positive
      # If a marker is not significant, but there is no adjacent QTL, declare a true negative
      # If a marker is signficant, but there is an adjacent QTL, declare a true positive
      # If a marker is not significant, but there is an adjacent marker, declare a false negative
      
      # Split by chromosome
      mar_assoc_qxe <- loci_sig_qxe %>% 
        map(~split(., .$chr) %>% 
              map(function(loci_map) {
                loci_map <- as.data.frame(loci_map) %>%
                  arrange(pos)
                
                assoc_results <- vector("list", nrow(loci_map))
              for (i in seq_along(assoc_results)) {
                if (loci_map[i,"type"] == "qtl") {
                  assoc_results[i] <- NA
                  
                } else if (is.na(loci_map[i,"p_value"])) {
                  
                  assoc_results[i] <- NA
                  
                } else {
                  
                  # IS there a QTL in the interval?
                  is_qtl <- any(loci_map[c(i - 1, i + 1),"type"] == "qtl", na.rm = T)
                  # Is there qxe effect
                  is_qxe <- rowSums(loci_map[c(i - 1, i + 1), c("env1", "env2")] != 0) > 0
                  # Sig level
                  sig_levels <- loci_map[i,c("fdr0.01_sig", "fdr0.05_sig", "fdr0.1_sig")]
                  
                  if (any(is_qtl & is_qxe, na.rm = T)) {
                    assoc_results[[i]] <- ifelse(sig_levels, "true_positive", "false_negative")
                    
                  } else {
                    assoc_results[[i]] <- ifelse(sig_levels, "false_positive", "true_negative")
                    
                  }
                  
                }
                
              }
                  
              cbind(select(loci_map, -contains("fdr")), 
                    `dimnames<-`(x = do.call("rbind", assoc_results), list(NULL, str_c("fdr", sig_cutoff))))
            }) %>%
            do.call("rbind", .) )
      
      
      
      ## Determine the number of true QTL detected
      # If a QTL is flanked by at least one significant marker, the QTL is discovered
      # Split by chromosome
      qtl_discover_main <- loci_sig_main %>% 
        map(~split(., .$chr) %>% 
              map(function(loci_map) {
              
              # Pull out markers and marker interval
              mar_map <- subset(x = loci_map, type == "marker")
              mar_interval <- mar_map$pos
              # Marker significance
              mar_sig <- mar_map[,c("fdr0.01_sig", "fdr0.05_sig", "fdr0.1_sig")]
              
              # Determine if QTL were discovered
              loci_map_qtl <- subset(x = loci_map, type == "qtl") %>% 
                select(-contains("fdr"))  %>% 
                mutate(first_mar = findInterval(x = pos, vec = mar_interval), 
                       sec_mar = first_mar + 1)
              
              # If no QTL are on the chromosome, skip
              if (nrow(loci_map_qtl) > 0) {
                discovered <- t(apply(loci_map_qtl[,c("first_mar", "sec_mar")], 
                                      MARGIN = 1, FUN = function(qtl) 
                                        colMeans(mar_sig[c(qtl[1], qtl[2]),]) > 0)) 
                
              } else {
                discovered <- NULL
                
              }
                
              cbind(loci_map_qtl, discovered) %>% select(marker:env2, contains("fdr"))  }) %>%
            bind_rows() )
      
      
      
      qtl_discover_qxe <- loci_sig_qxe %>%
        map(~split(., .$chr) %>% 
              map(function(loci_map) {
              
              # Pull out markers and marker interval
              mar_map <- subset(x = loci_map, type == "marker")
              mar_interval <- mar_map$pos
              # Marker significance
              mar_sig <- mar_map[,c("fdr0.01_sig", "fdr0.05_sig", "fdr0.1_sig")]
              
              # Determine if QTL were discovered
              loci_map_qtl <- loci_map %>%
                subset(type == "qtl" & (env1 != 0 | env2 != 0)) %>%
                select(-contains("fdr"))  %>% 
                mutate(first_mar = findInterval(x = pos, vec = mar_interval), 
                       sec_mar = first_mar + 1)
              
              # If no QTL are on the chromosome, skip
              if (nrow(loci_map_qtl) > 0) {
                discovered <- t(apply(loci_map_qtl[,c("first_mar", "sec_mar")], 
                                      MARGIN = 1, FUN = function(qtl) 
                                        colMeans(mar_sig[c(qtl[1], qtl[2]),]) > 0))
                
              } else {
                discovered <- NULL
                
              }
              
              cbind(loci_map_qtl, discovered) %>% select(marker:env2, contains("fdr"))  }) %>%
            bind_rows() )
      
      
      
      
      # Separate true/false and negative/positive
      mar_assoc_main1 <- mar_assoc_main %>% 
        bind_rows() %>% 
        as_data_frame() %>% 
        gather(sig_level, result, -marker:-type) %>%
        filter(type == "marker")
      
      mar_assoc_qxe1 <- mar_assoc_qxe %>% 
        bind_rows() %>% 
        as_data_frame() %>% 
        gather(sig_level, result, -marker:-type) %>%
        filter(type == "marker")
      
      assoc_count <- bind_rows(mar_assoc_main1, mar_assoc_qxe1) %>%
        filter(!is.na(result)) %>%
        separate(col = result, into = c("true_false", "neg_pos"), sep = "_") %>% 
        group_by(model, sig_level, term, true_false, neg_pos) %>% 
        summarize(n = n()) %>%
        ungroup() %>%
        complete(model, sig_level, term, true_false, neg_pos)
      
      
      
      # Count the number of QTL that were discovered
      qtl_discover_main1 <- list(qtl_discover_main, models) %>% 
        pmap_df(~mutate(.x, model = .y, term = "main_effect")) %>% 
        as_data_frame() %>% 
        gather(sig_level, discovered, -marker:-env2, -model, -term) %>%
        mutate(discovered = as.logical(discovered))
      
      qtl_discover_qxe1 <- list(qtl_discover_qxe, models) %>% 
        pmap_df(~mutate(.x, model = .y, term = "qxe")) %>% 
        as_data_frame() %>%
        gather(sig_level, discovered, -marker:-env2, -model, -term) %>%
        mutate(discovered = as.logical(discovered))
      
      
      
      # Count the number of QTL that were discovered
      discovered_count <- bind_rows(qtl_discover_main1, qtl_discover_qxe1) %>%
        group_by(model, term, sig_level) %>% 
        summarize(n_discovered = sum(discovered, na.rm = T))
      
      # Return the association and discovery results
      list(assoc_count = assoc_count, discovered_count = discovered_count) }, .to = "data")
    
  }, mc.cores = n_cores)
    


# Save the results
## Determine if a file is already there
save_file <- file.path(result_dir, "S2MET_mapping_gwas_qxe_simulation.RData")
present <- file.exists(save_file)

# Use while loop to edit the file until it does not exist
i = 1
while(present) {
  save_file <- file.path(result_dir, str_c("S2MET_mapping_gwas_qxe_simulation", i, ".RData"))
  present <- file.exists(save_file)
}

save("sim_results", file = save_file)

