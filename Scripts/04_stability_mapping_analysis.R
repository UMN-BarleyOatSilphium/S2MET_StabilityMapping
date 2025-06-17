## Analysis of GWAS results
## 
## This script will outline the analysis of GWAS mapping procedures using the S2MET data to
## identify marker-trait associations for the mean per se of traits and phenotypic stability
## This script also includes code for generating figures related to this analysis
## 
## 


# Load other packages
library(LDheatmap)
library(GenomicFeatures)
library(GenomicRanges)

# Project and other directories
source(here::here("startup.R"))


# Significance threshold
fdr_level <- 0.10

# Path to the barley GFF file
barley_gff <- file.path(geno_dir, "Hordeum_vulgare.IBSC_v2.49.gff3.gz")

# According to Hamblin 2010, LD decays to around r^2 = 0.2 at about 3 cM.
# With this in mind, calculate the average Mb distance of 3 cM.
# Using the smoothed recombination data from T. Kono, this works
# out to about 4.6 Mb / 3 cM.

# Set a window outside of each significant SNP to search for genes
sig_snp_window <- 4.6 * 1e6







# Load data ----------------------------------------------------

# Load the mapping results
load(file.path(result_dir, "stability_gwas_out.RData"))
# Load the mapping with resampling results
load(file.path(result_dir, "stability_samples_gwas_out.RData"))




# trait, parameter, PC combinations
trait_param_pc <- univariate_gwas_farm_bestModel_qq %>%
  select(trait, parameter, nPC = bestModel_qq_nPC)

# Load the GFF file
barley_txdb <- makeTxDbFromGFF(file = barley_gff, format = "gff3")
barley_txdb2 <- ape::read.gff(file = barley_gff)

# Subset genes and regular chromosomes
barley_txdb2_genes <- barley_txdb2 %>% 
  filter(seqid %in% paste0("chr", 1:7, "H"), type == "gene") %>%
  # Edit the attributes
  mutate(attributes1 = str_split(attributes, ";"),
         attributes1 = map(attributes1, ~{
           x_split <- str_split(.x, "=")
           as_tibble(setNames(object = map(x_split, 2), map_chr(x_split, 1)))
         }))

barley_txdb2_genes_grange <- barley_txdb2_genes %>%
  as_tibble() %>% 
  select(-attributes) %>% 
  unnest(attributes1) %>% 
  GRanges()


## Read in the T3 QTL GFF information
# Find the files
t3_qtl_gff_files <- list.files(path = geno_dir, pattern = "T3", full.names = TRUE)

# Read in the files
t3_qtl_gff_txdb <- map_df(t3_qtl_gff_files, ~ape::read.gff(file = .x)) %>%
  as_tibble() %>%
  # Listify the attributes
  mutate(attributes = str_split(attributes, ";"),
         t3_trait = map_chr(attributes, ~str_split(str_subset(string = .x, pattern = "_Trait"), "=")[[1]][2]),
         trial = map(attributes, ~str_remove(string = str_extract(string = str_subset(string = .x, pattern = "Phenotype_trial"), 
                                                                      pattern = 'trial_code%3D[A-Za-z0-9_]*'), pattern = "trial_code%3D")),
         experiment = map_lgl(attributes, ~any(str_detect(.x, "Phenotype_exp"))),
         phenotyping_source = map2_chr(trial, experiment, ~ifelse(.y, "experiment_set", .x)))



t3_qtl_granges <- t3_qtl_gff_txdb %>%
  select(-trial, -experiment) %>%
  GRanges(.)







# Find the significant markers for each trait -----------------------------

univariate_gwas_sigmar <- inner_join(trait_param_pc, univariate_gwas_farm) %>%
  filter(trait %in% traits) %>%
  mutate(fdr_p_thresh = ifelse(is.na(fdr_p_thresh), -Inf, fdr_p_thresh),
         sigmar = map2(scores, fdr_p_thresh, ~filter(.x, P.value <= .y))) %>%
  select(-betapc)

# Number of significant markers per trait and parameter
univariate_gwas_sigmar %>%
  mutate(nSigmar = map_dbl(sigmar, nrow)) %>%
  select(trait, parameter, nSigmar)

# Create genomic ranges for the significant markers
univariate_gwas_sigmar_grange_info <- univariate_gwas_sigmar %>%
  select(trait, parameter, sigmar) %>%
  unnest(sigmar) %>%
  mutate(chrom = paste0("chr", chrom, "H"))

# Number of total significant markers
nrow(univariate_gwas_sigmar_grange_info)


univariate_gwas_sigmar_grange <- GRanges(
  seqnames = univariate_gwas_sigmar_grange_info$chrom, 
  ranges = IRanges(start = univariate_gwas_sigmar_grange_info$pos, end = univariate_gwas_sigmar_grange_info$pos),
  marker = univariate_gwas_sigmar_grange_info$marker,
  trait = univariate_gwas_sigmar_grange_info$trait,
  parameter = univariate_gwas_sigmar_grange_info$parameter,
  p_value = univariate_gwas_sigmar_grange_info$P.value,
  maf = univariate_gwas_sigmar_grange_info$maf,
  effect = univariate_gwas_sigmar_grange_info$effect
)


# Find the overlap between the sigmar regions and barley genes

# Merge by overlaps
univariate_gwas_sigmar_grange_geneOverlap <- mergeByOverlaps(query = univariate_gwas_sigmar_grange, subject = barley_txdb2_genes_grange, 
                                                             ignore.strand = TRUE, maxgap = sig_snp_window)


# Filter for named genes
univariate_gwas_sigmar_grange_geneOverlap %>%
  subset(!is.na(Name))


# Find overlap with trial-specific QTL
univariate_gwas_overlap_perTrialT3QTL <- mergeByOverlaps(query = univariate_gwas_sigmar_grange, subject = subset(t3_qtl_granges, phenotyping_source != "experiment_set"), 
                                                         ignore.strand = TRUE, maxgap = sig_snp_window)

# Find overlap with experiment-wide QTL
univariate_gwas_overlap_experimentSetT3QTL <- mergeByOverlaps(query = univariate_gwas_sigmar_grange, subject = subset(t3_qtl_granges, phenotyping_source == "experiment_set"), 
                                                         ignore.strand = TRUE, maxgap = sig_snp_window)






# Calculate discovery rate in the sampling scheme -------------------------

univariate_gwas_sigmar_grange_info

# Create a data.frame of all of the significant markers discovered in each 
# sampling rep
# 
# Combine the sample results with the sigmar
univariate_gwas_sample_sigmar <- univariate_gwas_sample_farm %>%
  filter(trait %in% traits) %>%
  # Remove the null results
  filter(!sapply(gwas_out, is.null)) %>%
  # Iterate over rows
  group_by(trait, parameter, pEnv, rep) %>%
  do(discovered_sig_markers = {
    row <- .
    
    # FDR p value threshold
    fdr_p <- ifelse(is.na(row$gwas_out[[1]]$fdr_p_thresh), -Inf, row$gwas_out[[1]]$fdr_p_thresh)
    
    # List the significant markers from the sample
    # return this
    row$gwas_out[[1]]$GWAS %>%
      filter(P.value <= fdr_p)
    
  }) %>% ungroup()


# Calculate the rate of discovery for all other markers (even those not discovered in the large analysis)

# First create a list of all markers discovered in the samples
univariate_gwas_samples_all_discovered_markers <- univariate_gwas_sample_sigmar %>%
  filter(trait %in% traits) %>%
  unnest(discovered_sig_markers) %>%
  distinct(trait, parameter, pEnv, marker) %>% 
  # Add significant marker from the original analysis
  left_join(., select(univariate_gwas_sigmar, trait, parameter, sigmar)) %>%
  mutate(orig_signif_marker = map2_lgl(marker, sigmar, ~.x %in% .y$marker)) %>%
  select(-sigmar)

# Determine the discovery rate of these marker
univariate_gwas_sample_discovery_rate <- univariate_gwas_samples_all_discovered_markers %>%
  group_by(trait, parameter, pEnv) %>%
  nest() %>%
  ungroup() %>%
  left_join(., univariate_gwas_sample_sigmar) %>%
  mutate(data = map2(data, discovered_sig_markers, ~mutate(.x, discovered = as.numeric(marker %in% .y$marker)))) %>%
  select(-discovered_sig_markers) %>% 
  unnest(data) %>% 
  group_by(trait, parameter, pEnv, marker, orig_signif_marker) %>% 
  summarize(discovery_rate = mean(discovered), .groups = "drop")
  


# Plot the rate of discovery
univariate_gwas_sample_discovery_rate %>%
  ggplot(aes(x = pEnv, y = discovery_rate, color = orig_signif_marker, group = marker)) +
  geom_line(linewidth = 0.8) +
  facet_grid(trait ~ parameter)



# Calculate the average marker discovery rate
univariate_gwas_sample_discovery_rate_avg <- univariate_gwas_sample_discovery_rate %>%
  filter(orig_signif_marker) %>%
  group_by(trait, parameter, pEnv) %>%
  summarize(avg_discovery_rate = mean(discovery_rate), .groups = "drop")
  


# Save the results
save("univariate_gwas_sample_discovery_rate", "univariate_gwas_sample_discovery_rate_avg", 
     "univariate_gwas_sigmar_grange" , "t3_qtl_granges",
     "univariate_gwas_overlap_perTrialT3QTL", "univariate_gwas_overlap_experimentSetT3QTL",
     "univariate_gwas_sigmar_grange_geneOverlap", 
     file = file.path(result_dir, "univariate_gwas_analysis.RData"))





