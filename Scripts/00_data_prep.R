## S2MET StabilityMapping
##
## Data preparation
##
##


# Load packages and paths
source("startup.R")

library(snps) # Personal package for SNP data curation


## User parameters
max_LD_r2 <- 0.99999 # Max LD of any one pair of markers
max_LD_r2_2 <- 0.95
min_mac <- 10 # Minimum minor allele count
max_indiv_miss <- 0.7
max_snp_miss <- 0.7



# Load data ---------------------------------------------------------------


# Load the phenotypic data
load(file.path(pheno_dir, "S2_tidy_BLUE.RData"))

# Load the trial metadata
trial_info <- read_csv(file = file.path(data_dir, "trial_metadata.csv")) %>%
  filter(project2 == "S2MET") %>%
  ## Replace Ithaca1 and Ithaca2 with Ithaca
  mutate(location = str_replace_all(location, "Ithaca1|Ithaca2", "Ithaca"))

n_distinct(trial_info$environment)


# Load the genotypic data
load(file.path(geno_dir, "S2_genos_mat.RData"))

# Load BOPA snps
bopa_snps <- read_table(file = file.path(data_dir, "BCAP_Data/genotype.hmp.txt"))
# Convert to marker matrix
bopa_snp_info <- bopa_snps %>%
  select(marker = 1, allele = alleles, chrom, pos) %>%
  filter(! chrom %in% c("Un", "UNK")) %>%
  mutate(chrom = parse_number(chrom))
bopa_snp_mat <- bopa_snps %>%
  select(marker = 1, names(.), -alleles:-pos) %>%
  as.data.frame() %>%
  column_to_rownames("marker") %>%
  as.matrix() %>%
  t()



# Manage entry information ------------------------------------------------

# Load an entry file
entry_list <- read_excel(file.path(data_dir, "project_entries.xlsx"))


# Grab the entry names that are not checks
tp <- subset(entry_list, Class == "S2TP", Line, drop = TRUE)
tp <- setdiff(tp, c("07MT-10")) # Remove hullness line
# Program names for the tp
tp_program_name <- subset(entry_list, Line %in% tp) %>% 
  {setNames(object = .$Program_Name, nm = .$Line)}

# Grab the vp
vp <- subset(entry_list, Class == "S2C1R", Line, drop = TRUE)

# Vectors of tp and vp that are in the genotype matrix
tp_geno <- intersect(tp, row.names(s2_discrete_mat))
vp_geno <- intersect(vp, row.names(s2_discrete_mat))


# Manage genotype data  ---------------------------------------------------

# Construct the K matrix using all SNPs and all individuals
# Reduce to the FP and OP
Kall <- A.mat(X = s2_imputed_mat, min.MAF = 0, max.missing = 1)

Ksnp <- Kall[c(tp_geno, vp_geno), c(tp_geno, vp_geno)]

## Add missing entries as unrelated
Kmiss_tp <- diag(x = mean(diag(Ksnp[tp_geno, tp_geno])), nrow = length(tp) - length(tp_geno))
dimnames(Kmiss_tp) <- replicate(2, setdiff(tp, tp_geno), simplify = FALSE)
Kmiss_vp <- diag(x = mean(diag(Ksnp[vp_geno, vp_geno])), nrow = length(vp) - length(vp_geno))
dimnames(Kmiss_vp) <- replicate(2, setdiff(vp, vp_geno), simplify = FALSE)
# Combine these
Kmiss <- as.matrix(Matrix::bdiag(Kmiss_tp, Kmiss_vp))

# Combine with the Ksnp matrix
K <- as.matrix(Matrix::bdiag(Ksnp, Kmiss))
dimnames(K) <- list(c(row.names(Ksnp), row.names(Kmiss_tp), row.names(Kmiss_vp)), c(row.names(Ksnp), row.names(Kmiss_tp), row.names(Kmiss_vp)))



## Marker filtered for GWAS

# Combine GBS and BOPA SNPs; filter for the TP
snp_info_combined <- bind_rows(rename(snp_info, marker = 1), bopa_snp_info) %>%
  arrange(chrom, pos)
  
# Rename BOPA genotypes to the TCAP name
tp_program_geno <- subset(tp_program_name, tp_program_name %in% row.names(bopa_snp_mat))
bopa_snp_mat_tp <- bopa_snp_mat[tp_program_geno,, drop = FALSE]
row.names(bopa_snp_mat_tp) <- names(tp_program_geno)

# Expand each matrix to include missing data
tp_geno_all <- union(tp_geno, names(tp_program_geno))
  
gbs_matrix_expand <- matrix(data = NA, nrow = length(setdiff(tp_geno_all, tp_geno)), ncol = ncol(s2_imputed_mat),
                            dimnames = list(setdiff(tp_geno_all, tp_geno), colnames(s2_imputed_mat)))
s2_imputed_mat1 <- rbind(s2_imputed_mat[c(tp_geno), , drop = FALSE], gbs_matrix_expand)

bopa_matrix_expand <- matrix(data = NA, nrow = length(setdiff(tp_geno_all, names(tp_program_geno))), ncol = ncol(bopa_snp_mat_tp),
                             dimnames = list(setdiff(tp_geno_all, names(tp_program_geno)), colnames(bopa_snp_mat_tp)))
bopa_snp_mat_tp1 <- rbind(bopa_snp_mat_tp, bopa_matrix_expand)

# Combine
geno_mat_tp <- cbind(s2_imputed_mat1, bopa_snp_mat_tp1)[,snp_info_combined$marker, drop = FALSE]

# Calculate minor allele frequency; plot
maf1 <- calc_maf2(x = geno_mat_tp); hist(maf1, main = "MAF - prefilter")

# For any overlapping marker (in the exact same physical positions), retain the marker
# with the higher MAF
overlap_marker_info <- snp_info_combined %>%
  group_by(chrom, pos) %>% 
  filter(n() > 1) %>%
  left_join(., tibble(marker = names(maf1), maf = maf1)) %>%
  ungroup()

# Pick the marker with higher MAF
overlap_marker_retained <- overlap_marker_info %>%
  arrange(chrom, pos, desc(maf)) %>%
  group_by(chrom, pos) %>%
  slice(1) %>%
  ungroup()


geno_mat_tp1 <- geno_mat_tp[, setdiff(colnames(geno_mat_tp), setdiff(overlap_marker_info$marker, overlap_marker_retained$marker))]


# Filter on MAF - this will be used to calculate LD
geno_mat_filter2 <- filter_snps2(x = geno_mat_tp1, r2.max = max_LD_r2, maf.min = min_mac / nrow(geno_mat_tp1),
                                 indiv.miss.max = 0.8, snp.miss.max = max_snp_miss)
geno_mat_LD <- geno_mat_filter2

# Dimensions
dim(geno_mat_LD)

# Now filter on MAF AND a lower LD threshold - these will be used for mapping
geno_mat_filter3 <- filter_snps2(x = geno_mat_filter2, r2.max = max_LD_r2_2, maf.min = min_mac / nrow(geno_mat_filter2),
                                 indiv.miss.max = 0.8, snp.miss.max = max_snp_miss)
geno_mat_filter3

# Impute
geno_mat_mapping_out <- A.mat(X = geno_mat_filter3, min.MAF = 0, max.missing = 1, return.imputed = TRUE, impute.method = "EM")
geno_mat_mapping <- geno_mat_mapping_out$imputed

# Dimensions
dim(geno_mat_mapping)

# Recalculate and visualize MAF
maf2 <- calc_maf2(x = geno_mat_mapping); hist(maf2, main = "MAF - postfilter")

# Copy snp info
snp_info_all <- snp_info_combined

# filter the snp info
snp_info <- snp_info_combined %>%
  filter(marker %in% colnames(geno_mat_mapping))


## Remove environments with low heritability
## This will be the df for correcting for heritability when calculating prediction accuracy
env_trait_herit <- s2_metadata %>% 
  filter(trait %in% traits, trial %in% trial_info$trial, heritability >= 0.10) %>%
  select(trial, trait, heritability)


# Filter BLUEs
pheno_dat <- s2_tidy_BLUE %>%
  filter(trial %in% trial_info$trial,
         line_name %in% c(tp, vp),
         trait %in% traits) %>%
  # Filter out environments with low heritability
  inner_join(., select(env_trait_herit, -heritability)) %>%
  # Add full location name
  select(-location) %>%
  left_join(., distinct(trial_info, environment, location)) %>%
  # Remove environments deemed failures (i.e. HNY16 for grain yield)
  filter(!(environment == "HNY16" & trait == "GrainYield"),
         !(environment == "EON17" & trait == "HeadingDate"),
         !(environment == "KNY16" & trait == "TestWeight")) %>%
  # Rename and reorder
  select(trial, environment, location, year, trait, line_name, value, std_error = std.error)


# Save everything
save("pheno_dat", "trial_info", "snp_info", "snp_info_all", "geno_mat_mapping", 
     "geno_mat_LD", "K", "tp", "vp", "tp_geno", "vp_geno",
     file = file.path(data_dir, "project_data.RData"))

