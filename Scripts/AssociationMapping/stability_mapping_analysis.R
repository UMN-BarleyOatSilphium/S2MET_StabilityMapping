## Analysis of GWAS results
## 
## This script will outline the analysis of GWAS mapping procedures using the S2MET data to
## identify marker-trait associations for the mean per se of traits and phenotypic stability
## This script also includes code for generating figures related to this analysis
## 
## Author: Jeff Neyhart
## Last updated: June 27, 2018
## 


###########################################################3
## Run the following for all analyses


library(qvalue)
library(cowplot)
library(broom)
library(LDheatmap)
library(GenomicRanges)


# Run the source script
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))


# Rename the marker matrix
M <- S2TP_imputed_multi_genos_mat

# Get the marker information
snp_info <- S2TP_imputed_multi_genos_hmp %>%
  select(marker = rs, chrom, pos, cM_pos) %>%
  # Correct the cM position for BOPA snps
  mutate(cM_pos = if_else(str_detect(marker, "^S"), cM_pos, cM_pos / 1000))

# Subset the entry list for the tp
tp_entry_list <- entry_list %>%
  filter(Line %in% tp_geno) %>%
  select(line_name = Line, program = Program)

# Significance threshold for GWAS
alpha <- 0.05

# What model will we use for mapping?
model_use <- "QG"


# Define a window around the significant SNPs
window <- 5e6 # 5 Mbp


#######################################################






### Population Structure

# First examine population structure within the germplasm

# Calculate the kinship matrix
K <- A.mat(X = M, min.MAF = 0, max.missing = 1)

## Visualize pop structure
K_prcomp <- prcomp(K)

# Tidy, calculate lambda, then add program information
K_prcomp_df <- tidy(K_prcomp) %>%
  mutate(PC = str_c("PC", PC)) %>%
  dplyr::rename(line_name = row) %>%
  left_join(entry_list, by = c("line_name" = "Line"))

# Extract lambda and calculate variance explained
var_exp <- K_prcomp$sdev %>% 
  set_names(str_c("PC", seq_along(.))) %>%
  {. / sum(.)}

## Combine data.frame to plot
df_PC1_PC2 <- K_prcomp_df %>% 
  filter(PC %in% c("PC1", "PC2")) %>% 
  spread(PC, value) %>%
  mutate(x = "PC1", y = "PC2", var_expx = var_exp["PC1"], var_expy = var_exp["PC2"]) %>%
  dplyr::rename(xvalue = PC1, yvalue = PC2)

df_PC1_PC3 <- K_prcomp_df %>% 
  filter(PC %in% c("PC1", "PC3")) %>% 
  spread(PC, value) %>%
  mutate(x = "PC1", y = "PC3", var_expx = var_exp["PC1"], var_expy = var_exp["PC3"]) %>%
  dplyr::rename(xvalue = PC1, yvalue = PC3)

df_PC2_PC3 <- K_prcomp_df %>% 
  filter(PC %in% c("PC2", "PC3")) %>% 
  spread(PC, value) %>%
  mutate(x = "PC2", y = "PC3", var_expx = var_exp["PC2"], var_expy = var_exp["PC3"]) %>%
  dplyr::rename(xvalue = PC2, yvalue = PC3)

## Combine the data frames and combine the variance explained with
## the PC
df_combined <- bind_rows(df_PC1_PC2, df_PC1_PC3, df_PC2_PC3) %>%
  mutate(x1 = str_c(x, " (", 100 * round(var_expx, 3), "%)"),
         y1 = str_c(y, " (", 100 * round(var_expy, 3), "%)"))


colors_use <- set_names(umn_palette(2)[3:7], unique(df_combined$Program))

## Plot
g_pop_str <- df_combined %>% 
  ggplot(aes(x = xvalue, y = yvalue, col = Program)) + 
  geom_point(size = 0.5) + 
  facet_grid(y1 ~ x1, switch = "both") +
  scale_color_manual(name = "Breeding\nProgram", values = colors_use) +
  theme_pnas() +
  theme(legend.position = c(0.75, 0.75), axis.title.x = element_blank(),
        legend.key.height = unit(0.75, units = "lines"))

# Save this
ggsave(filename = "population_structure.jpg", plot = g_pop_str, path = fig_dir,
       width = 8.7, height = 8, units = "cm", dpi = 1000)





## How is population structure correlated with the traits?

# Load the genotype means and FW regression results
load(file.path(result_dir, "pheno_mean_fw_results.RData"))

pheno_mean_fw <- pheno_mean_fw_tpvp %>%
  filter(line_name %in% tp)


# Transform the delta statistic to log-delta
pheno_mean_fw_trans <- pheno_mean_fw %>% 
  distinct(trait, line_name, g, b, delta) %>% 
  mutate(log_delta = log(delta)) %>%
  select(-delta)


# Combine the trait and PCs, duplicating the BLUEs to make summarizing smoother
trait_pop_str <- K_prcomp_df %>% 
  filter(PC %in% c("PC1", "PC2", "PC3")) %>% 
  select(line_name, program = Program, PC, eigenvalue = value) %>%
  left_join(., pheno_mean_fw_trans) %>%
  gather(measure, value, g:log_delta)


## Correlate the eigenvalues with the trait mean and stability, and calculate regression coefficients
## Use a bootstrap to estimate confidence interval
trait_pop_str_corr <- trait_pop_str %>% 
  group_by(trait, PC, measure) %>% 
  do(neyhart::bootstrap(x = .$value, y = .$eigenvalue, fun = "cor", boot.reps = 1000)) %>%
  # Which ones are significant
  mutate(significant = !between(0, ci_lower, ci_upper),
         sig_ann = ifelse(significant, "*", ""),
         annotation = str_c("r = ", formatC(base, digits = 3, format = "f"), sig_ann))

## Parametric correlation test
trait_pop_str_corr_param <- trait_pop_str %>% 
  group_by(trait, PC, measure) %>% 
  do(test = cor.test(x = .$value, y = .$eigenvalue)) %>%
  ungroup() %>%
  mutate(base = map_dbl(test, "estimate"),
         pvalue = map_dbl(test, "p.value"),
         significant = pvalue <= alpha,
         sig_ann = ifelse(significant, "*", ""),
         annotation = str_c("r = ", formatC(base, digits = 3, format = "f"), sig_ann))


trait_pop_str_annotation <- trait_pop_str_corr %>%
  ungroup() %>%
  distinct(trait, measure, PC, annotation) %>%
  mutate(measure = str_replace_all(measure, coef_replace))

# Add the correlations back to the data
trait_pop_str1 <- left_join(trait_pop_str, trait_pop_str_corr)

## Common plot modifier
g_mod <- list(
  geom_point(),
  geom_smooth(method = "lm", se = FALSE, color = "black"),
  geom_text(aes(x = Inf, y = -Inf, label = annotation), hjust = 1.5, vjust = 0.25, color = "black", size = 6),
  facet_grid(trait ~ PC, scales = "free"),
  scale_color_discrete(name = "Breeding\nProgram"),
  xlab("Eigenvalue"),
  theme_bw())


# Plot pop str vs mean
g_PC_v_mean <- trait_pop_str1 %>% 
  filter(measure == "g") %>%
  ggplot(aes(x = eigenvalue, y = value, col = program)) + 
  ylab("Genotypic Value") +
  g_mod

# Plot pop str vs linear stability
g_PC_v_b <- trait_pop_str1 %>% 
  filter(measure == "b") %>%
  ggplot(aes(x = eigenvalue, y = value, col = program)) + 
  ylab("Linear Stability") +
  g_mod

# Plot pop str vs nonlinear stability
g_PC_v_delta <- trait_pop_str1 %>%
  filter(measure == "log_delta") %>%
  ggplot(aes(x = eigenvalue, y = value, col = program)) + 
  ylab("Non-Linear Stability") +
  g_mod

# Save
ggsave(filename = "population_structure_versus_mean.jpg", plot = g_PC_v_mean,
       path = fig_dir, height = 4, width = 6, dpi = 1000)  

ggsave(filename = "population_structure_versus_linear_stability.jpg", plot = g_PC_v_b,
       path = fig_dir, height = 4, width = 6, dpi = 1000)  

ggsave(filename = "population_structure_versus_nonlinear_stability.jpg", plot = g_PC_v_delta,
       path = fig_dir, height = 4, width = 6, dpi = 1000)  


## Just PC1
g_PC1_v_traits <- trait_pop_str1 %>% 
  filter(PC == "PC1") %>% 
  mutate(measure = str_replace_all(measure, coef_replace)) %>%
  ggplot(aes(x = eigenvalue, y = value, col = program)) + 
  geom_point(size = 0.2) + 
  geom_text(data = subset(trait_pop_str_annotation, PC == "PC1"), aes(x = Inf, y = Inf, label = annotation), 
            inherit.aes = FALSE, hjust = 1, vjust = 1.5, size = 2) + 
  geom_smooth(method = "lm", se = FALSE, color = "black", lwd = 0.2) +
  scale_color_manual(values = colors_use, name = "Breeding\nProgram") +
  facet_wrap(trait ~ measure, scales = "free") +
  xlab("PC1") +
  theme_pnas() + 
  theme(axis.title.y = element_blank())

## Combine plots

g_pop_str_fig <- plot_grid(
  g_pop_str + ylab(""),
  g_PC1_v_traits + theme(legend.position = "none"),
  ncol = 2, labels = LETTERS[1:2], align = "hv", axis = "lr", hjust = 0.001
)

# Save
ggsave(filename = "population_structure_and_traits.jpg", plot = g_pop_str_fig, path = fig_dir,
       width = 17.8, height = 8, unit = "cm", dpi = 1000)










### Association Results - Genomewide Scan

## Load the results of the genomewide scan
load(file.path(result_dir, "pheno_fw_mean_gwas_results.RData"))

# Use the results from the TP+VP analysis
gwas_pheno_mean_fw_tidy_adj <- gwas_pheno_mean_fw_tpvp_tidy_adj
gwas_mlmm_final <- gwas_mlmm_final_tpvp


## Select the desired model
gwas_pheno_mean_fw_adj <- subset(gwas_pheno_mean_fw_tidy_adj, model == model_use)

gwas_mlmm_model_adj <- subset(gwas_mlmm_final, model == model_use) %>%
  left_join(., select(snp_info, marker:pos), c("term" = "marker")) %>%
  select(trait:model, marker = term, chrom:pos, beta, pvalue, qvalue, snp_r_squared = r_squared,
          all_snp_r_squared = all_snps)
  


#### Summary of significant associations

# Sumarize
gwas_pheno_mean_fw_sig %>% 
  group_by(trait, coef) %>% 
  summarize_at(vars(contains("SNP")), sum)

# trait       coef      GWAS_sig_SNP MLMM_sig_SNP
# 1 GrainYield  b                    0            0
# 2 GrainYield  g                    0            0
# 3 GrainYield  log_delta            0            0
# 4 HeadingDate b                   20            4
# 5 HeadingDate g                   10            2
# 6 HeadingDate log_delta            3            3
# 7 PlantHeight b                   49            2
# 8 PlantHeight g                    3            2
# 9 PlantHeight log_delta            0            0

## GWAS: 85
## MLMM: 13

# Look at the significant markers and assess the effect size, minor allele frequency, etc.

# Calculate the allele frequency of the 1 allele
af1 <- data.frame(af = colMeans(M + 1) / 2 ) %>% 
  rownames_to_column("marker") %>% 
  mutate(maf = pmin(af, 1 - af))

# Combine the minor allele frequency information
gwas_mlmm_marker_info <- gwas_mlmm_model_adj %>% 
  ungroup() %>%
  left_join(., af1) %>%
  select(-af)

# Do the same for the genomewide scan
gwas_sig_marker_info <- gwas_pheno_mean_fw_adj %>%
  filter(qvalue <= alpha) %>% 
  select(trait, coef, model, names(.), -color) %>%
  left_join(., af1) %>%
  select(-af)



# Subset the marker matrix for these SNPs and count the number of lines
# from each program with each allele
gwas_mlmm_marker_prop <- gwas_mlmm_marker_info %>%
  group_by(marker) %>%
  do({
    df <- .
    
    # Extract the data on the marker
    marker_geno <- M[,df$marker, drop = FALSE] %>% 
      as.data.frame() %>% rownames_to_column("line_name") %>% 
      mutate(program = str_extract(line_name, "[A-Z][A-Z0-9]")) %>%
      select(line_name, program, marker = df$marker)
    
    # Round the marker genotype
    marker_geno_round <- marker_geno %>% 
      mutate(allele = case_when(marker <= -0.5 ~ -1,
                                marker >= 0.5 ~ 1,
                                TRUE ~ 0),
             allele = allele + 1)
    
    # Calculate the % of lines from each program that have the 1 allele
    marker_geno_prop <- marker_geno_round %>% 
      group_by(program) %>% 
      summarize(prop_1_allele = mean(allele) / 2) %>% 
      spread(program, prop_1_allele) 
    
    # Return the df
    bind_cols(df, marker_geno_prop)
    
  }) %>% arrange(trait, coef, chrom, pos)



## Pretty manhattan plots
# Add the position of the preceding chromosome to the position of the next
chrom_pos_cumsum <- snp_info %>% 
  group_by(chrom) %>% 
  summarize(chrom_pos = max(pos)) %>% 
  mutate(max_pos = c(0, head(chrom_pos, -1)), 
         next_chrom_start_pos = cumsum(max_pos),
         new_chrom_end = chrom_pos + next_chrom_start_pos,
         chrom_label_pos =( next_chrom_start_pos + new_chrom_end) / 2)

snp_info_new_pos <- left_join(snp_info, select(chrom_pos_cumsum, -max_pos)) %>% 
  mutate(new_pos = pos + next_chrom_start_pos)

  
g_all_gwas <- gwas_pheno_mean_fw_adj %>% 
  left_join(snp_info_new_pos) %>%
  # filter(coef != "log_delta") %>%
  mutate(coef = str_replace_all(coef, coef_replace)) %>%
  ggplot(aes(x = new_pos, y = neg_log_p, color = coef, shape = coef)) + 
  # Add chromosome lines
  geom_segment(data = chrom_pos_cumsum, aes(x = new_chrom_end, xend = new_chrom_end, y = 0, yend = 10), inherit.aes = FALSE, color = "grey75") +
  geom_point() +
  geom_text(data = chrom_pos_cumsum, aes(x = chrom_label_pos, y = -1, label = chrom), vjust = 0, inherit.aes = FALSE) +
  scale_color_discrete(name = NULL) +
  scale_shape_discrete(name = NULL) +
  ylab(expression(-log[10](italic(p)))) + 
  ylim(c(0, 10)) + 
  facet_grid(trait ~ ., switch = "y", space = "free_x") + 
  theme_classic() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), 
        strip.background = element_rect(fill = "grey85", size = 0), legend.position = c(0.08, 0.9),
        legend.text = element_text(size = 6))


# Save
ggsave(filename = "gwas_manhattan_all_trait_QG.jpg", plot = g_all_gwas, path = fig_dir, 
       height = 6, width = 10, dpi = 1000)


  
  

### LD Heatmap

# # For each significant region, look at a LD heatmap of surrounding markers
# gwas_mlmm_LD <- gwas_mlmm_marker_prop %>%
#   ungroup() %>%
#   mutate(ld_heatmap = vector("list", nrow(.)))
# 
# for (i in seq(nrow(gwas_mlmm_LD))) {
#   
#   snp <- gwas_mlmm_LD[i,]
#   
#   # Subset these SNPs and calculate LD
#   surrounding_snps <- snp_info %>% 
#     filter(chrom %in% snp$chrom, between(pos, mean(snp$pos) - window, mean(snp$pos) + window))
#   surrounding_snp_LD <- LD(x = M[,surrounding_snps$marker, drop = FALSE], df = FALSE)
#   
#   # Get the qvalues for these SNPs
#   surrounding_snp_qvalue <- gwas_pheno_mean_fw_adj %>% 
#     filter(marker %in% surrounding_snps$marker, trait %in% snp$trait, coef %in% snp$coef)
#   
#   gwas_mlmm_LD$ld_heatmap[[i]] <- LDheatmap(surrounding_snp_LD, surrounding_snps$pos, flip = TRUE,
#                           SNP.name = snp$marker) %>%
#     LDheatmap.addScatterplot(LDheatmap = ., P = surrounding_snp_qvalue$neg_log_q)
#   
#   
# }

 

### Annotation
### 
 
# Gather significant markers and create a GRange object
gwas_mlmm_grange <- gwas_mlmm_marker_info %>%
  # Add twice the window to complement the other overlapping procedures
  mutate(start = pos - (window), end = pos + (window)) %>%
  makeGRangesFromDataFrame(., keep.extra.columns = TRUE)

gwas_sig_mar_grange <- gwas_sig_marker_info %>%
  # Add twice the window to complement the other overlapping procedures
  mutate(start = pos - (window), end = pos + (window)) %>%
  makeGRangesFromDataFrame(., keep.extra.columns = TRUE)



## Split the mlmm grange by trait and coefficient and find any markers that overlap
##
# Split by coefficient and create a GRangesList
gwas_mlmm_grange_split <- gwas_mlmm_grange %>%
  as.data.frame() %>%
  split(.$trait) %>%
  map(~split(., .$coef) %>% map(~makeGRangesFromDataFrame(df = ., keep.extra.columns = TRUE)) %>%
        GRangesList())

## Find overlaps, then remove the self overlaps
gwas_mlmm_overlaps_list <- gwas_mlmm_grange_split %>%
  map(~findOverlaps(query = ., subject = .)) %>%
  map(~subset(., queryHits != subjectHits)) %>%
  map(~as.data.frame(.))

# ## Explore the overlaps
# gwas_mlmm_overlaps <- list(gwas_mlmm_overlaps_list, as.list(gwas_mlmm_grange_split)) %>%
#   pmap(~apply(X = .x, MARGIN = 1, FUN = function(overlap) {
#     if (nrow(.x) == 0) {
#       return(NA)
#     } else {
#       mergeByOverlaps(query = as.list(.y)[[overlap[1]]], subject = as.list(.y)[[overlap[2]]])
#     }}) %>% as.list() %>% do.call("rbind", .)) %>%
#   do.call("rbind", .) %>%
#   subset(., , c(2, 3, 5, 16, 18)) %>%
#   as_data_frame()


## There were no MLMM markers for the mean that overlapped with stability



## Split the sig mar grange by trait and coefficient and find any markers that overlap
## 
# Split by coefficient and create a GRangesList
gwas_sig_mar_grange_split <- gwas_sig_mar_grange %>% 
  as.data.frame() %>%
  split(.$trait) %>%
  map(~split(., .$coef) %>% map(~makeGRangesFromDataFrame(df = ., keep.extra.columns = TRUE)) %>% 
        GRangesList())

## Find overlaps, then remove the self overlaps
gwas_sig_mar_overlaps_list <- gwas_sig_mar_grange_split %>%
  map(~findOverlaps(query = ., subject = .)) %>%
  map(~subset(., queryHits != subjectHits)) %>%
  map(~as.data.frame(.))

## Explore the overlaps
gwas_sig_mar_overlaps <- list(gwas_sig_mar_overlaps_list, as.list(gwas_sig_mar_grange_split)) %>%
  pmap(~apply(X = .x, MARGIN = 1, FUN = function(overlap) {
    if (nrow(.x) == 0) {
      return(NA)
    } else {
      mergeByOverlaps(query = as.list(.y)[[overlap[1]]], subject = as.list(.y)[[overlap[2]]])
    }}) %>% as.list() %>% do.call("rbind", .)) %>% 
  do.call("rbind", .) %>%
  subset(., , c(2, 3, 5, 13, 15)) %>%
  as_data_frame()










### Check the significant markers for overlap with loci previously associated with
### the mean per se of agronomic traits

# Read in the QTL metadata files that were downloaded from T3
qtl_meta <- map_df(list.files(data_dir, pattern = "qtl_meta", full.names = TRUE), read_csv)
# Read in association data from Pauli2014 and Wang2012
bcap_association <- read_csv(file.path(data_dir, "BCAP_association_qtl.csv")) %>%
  mutate(position = parse_number(position)) %>%
  select(trait, marker, chrom = chromosome, position, gene:feature, reference) %>%
  filter(trait %in% c("GrainYield", "HeadingDate", "PlantHeight"))

# Read in known genes
known_genes <- read_csv(file = file.path(data_dir, "known_genes.csv"))


# Rename the traits
qtl_meta_df <- qtl_meta %>%
  mutate(trait = case_when(trait == "grain yield" ~ "GrainYield",
                           trait == "heading date" ~ "HeadingDate",
                           trait == "plant height" ~ "PlantHeight"))

# Remove unmapped QTL and change the chromosome name
qtl_meta_use <- qtl_meta_df %>% 
  filter(!chromosome %in% c("chrUNK", "chrUn")) %>% 
  mutate(chrom = as.integer(parse_number(chromosome))) %>%
  select(trait, marker, chrom, position, gene, feature) %>%
  arrange(trait, chrom, position)

# Combine data
qtl_meta_use1 <- bind_rows(qtl_meta_use, bcap_association)



## Create GRanges for each marker data.frame
qtl_meta_grange <- qtl_meta_use1 %>%
  split(.$trait) %>%
  map(~makeGRangesFromDataFrame(., keep.extra.columns = TRUE, start.field = "position", 
                                end.field = "position")) 


## Find the overlap with the MLMM markers
common_traits <- intersect(names(gwas_mlmm_grange_split), names(qtl_meta_grange))

gwas_mlmm_qtl_meta_overlap <- list(gwas_mlmm_grange_split[common_traits], qtl_meta_grange[common_traits]) %>%
  pmap(., ~lapply(as.list(.x), function(mar) mergeByOverlaps(query = mar, subject = .y))) %>%
  map(~do.call("rbind", .)) %>% do.call("rbind", .) %>%
  subset(x = ., , c(trait, coef, marker, marker, gene, reference)) %>%
  as_data_frame()

## Find the overlap with all sig markers
common_traits <- intersect(names(gwas_sig_mar_grange_split), names(qtl_meta_grange))

gwas_sig_mar_qtl_meta_overlap <- list(gwas_sig_mar_grange_split[common_traits], qtl_meta_grange[common_traits]) %>%
  pmap(., ~lapply(as.list(.x), function(mar) mergeByOverlaps(query = mar, subject = .y))) %>%
  map(~do.call("rbind", .)) %>% do.call("rbind", .) %>%
  subset(x = ., , c(trait, coef, marker, marker, gene, reference)) %>%
  as_data_frame()



## Nest the marker overlaps and combine with the original data.frame for printing
gwas_mlmm_marker_info_ann <- gwas_mlmm_marker_info %>%
  left_join(., group_by(gwas_mlmm_qtl_meta_overlap, trait, coef, marker) %>% 
              nest(marker.1:reference, .key = "annotation")) %>%
  mutate(qtl_hits = map_dbl(annotation, ~ifelse(is.null(.), 0, nrow(.))))

## Count the significant markers with overlaps
(gwas_mlmm_marker_info_ann1 <- gwas_mlmm_marker_info_ann %>% 
  group_by(trait, coef) %>% 
  summarize(prop_qtl_hits = mean(!map_lgl(annotation, is.null)), 
            sum_qtl_hits = sum(qtl_hits)))

## These are the number of MLMM markers with overlap to known genes or QTL


# trait       coef      prop_qtl_hits sum_qtl_hits
# 1 HeadingDate b                 0.5              2
# 2 HeadingDate g                 1                3
# 3 HeadingDate log_delta         0.333            1
# 4 PlantHeight b                 1                5
# 5 PlantHeight g                 0.5              1

gwas_sig_marker_info_ann <- gwas_sig_marker_info %>%
  left_join(., group_by(gwas_mlmm_qtl_meta_overlap, trait, coef, marker) %>% 
              nest(marker.1:reference, .key = "annotation")) %>%
  mutate(qtl_hits = map_dbl(annotation, ~ifelse(is.null(.), 0, nrow(.))))
  
## Count the significant markers with overlaps
(gwas_sig_marker_info_ann1 <- gwas_sig_marker_info_ann %>% 
  group_by(trait, coef) %>% 
  summarize(prop_qtl_hits = mean(!map_lgl(annotation, is.null)), 
            sum_qtl_hits = sum(qtl_hits)))

# trait       coef      prop_qtl_hits sum_qtl_hits
# 1 HeadingDate b                0.1               2
# 2 HeadingDate g                0.2               3
# 3 HeadingDate log_delta        0.333             1
# 4 PlantHeight b                0.0408            5
# 5 PlantHeight g                0.333             1



## Output a table
gwas_mlmm_marker_toprint <- gwas_mlmm_marker_info_ann %>%
  select(Trait = trait, Coef = coef, Marker = marker, Chrom = chrom, Position = pos, 
         Beta = beta, pvalue = pvalue, R2 = snp_r_squared, total_R2 = all_snp_r_squared, MAF = maf, QTLHits = qtl_hits,
         annotation) %>% 
  arrange(Trait, Coef, Chrom, Position) %>% 
  mutate(Coef = str_replace_all(Coef, coef_replace),
         Reference = map(annotation, "reference") %>% map(~ifelse(is.na(.), "T3", .)) %>% 
           map_chr(~paste0(unique(.), collapse = ", "))) %>%
  select(-annotation) %>%
  # Round
  mutate_at(vars(Beta:MAF), funs(formatC(., digits = 3, format = "g")))

write_csv(x = gwas_mlmm_marker_toprint, path = file.path(fig_dir, "gwas_mlmm_significant_associations.csv"))


## Output a table
gwas_sig_marker_toprint <- gwas_sig_marker_info_ann %>%
  select(Trait = trait, Coef = coef, Marker = marker, Chrom = chrom, Position = pos, 
         qvalue = qvalue, MAF = maf, QTLHits = qtl_hits, annotation) %>% 
  arrange(Trait, Coef, Chrom, Position) %>% 
  mutate(Coef = str_replace_all(Coef, coef_replace),
         Reference = map(annotation, "reference") %>% map(~ifelse(is.na(.), "T3", .)) %>% 
           map_chr(~paste0(unique(.), collapse = ", "))) %>%
  select(-annotation)

write_csv(x = gwas_sig_marker_toprint, path = file.path(fig_dir, "gwas_significant_associations.csv"))


## Use the annotation in the GWAS

## Pretty manhattan plots
# Add the position of the preceding chromosome to the position of the next
chrom_pos_cumsum <- snp_info %>% 
  group_by(chrom) %>% 
  summarize(chrom_pos = max(pos)) %>% 
  mutate(max_pos = c(0, head(chrom_pos, -1)), 
         next_chrom_start_pos = cumsum(max_pos),
         new_chrom_end = chrom_pos + next_chrom_start_pos,
         chrom_label_pos =( next_chrom_start_pos + new_chrom_end) / 2)

snp_info_new_pos <- left_join(snp_info, select(chrom_pos_cumsum, -max_pos)) %>% 
  mutate(new_pos = pos + next_chrom_start_pos)

## Adjust annotation positions
qtl_meta_use1_newpos <- qtl_meta_use1 %>% 
  select(trait, chrom, position) %>% 
  left_join(chrom_pos_cumsum) %>% 
  mutate(new_qtl_pos = position + next_chrom_start_pos)

known_genes_newpos <- known_genes %>% 
  select(trait:end) %>% 
  left_join(chrom_pos_cumsum) %>% 
  mutate(new_gene_pos = ((start + end) / 2) + next_chrom_start_pos)

# Adjust the position of the label of some genes
known_genes_newpos1 <- known_genes_newpos %>% 
  mutate(label_pos = new_gene_pos, 
         label_pos = ifelse(gene %in% c("Ppd-H2", "denso", "eps4L", "Vrn-H1"), label_pos - 2e8, label_pos),
         label_pos = ifelse(gene %in% c("eps7S", "Vrn-H3"), label_pos + 2e8, label_pos)) %>%
  # Remove eps7S
  filter(!gene %in% c("eps7S", "eps2"))


## Reformat the MLMM results to highlight as large points in the manhattan plot
gwas_mlmm_marker_info_ann_toplot <- gwas_mlmm_marker_info_ann %>% 
  left_join(snp_info_new_pos) %>% 
  mutate(coef = str_replace_all(coef, coef_replace),
         neg_log_q = -log10(qvalue)) %>% 
  select(trait:pos, new_pos, neg_log_q)
  



# New colors
colors_use <- set_names(umn_palette(2)[3:5], coef_replace)


g_all_gwas_annotated <- gwas_pheno_mean_fw_adj %>% 
  left_join(snp_info_new_pos) %>%
  # filter(coef != "log_delta") %>%
  mutate(coef = str_replace_all(coef, coef_replace)) %>%
  # ggplot(aes(x = new_pos, y = neg_log_p, color = coef, shape = coef)) + 
  ggplot(aes(x = new_pos, y = neg_log_q, color = coef, shape = coef)) + 
  # Significance thresholds
  geom_hline(yintercept = -log10(alpha), lwd = 0.2, lty = 3) + 
  # Add chromosome lines
  # geom_segment(data = chrom_pos_cumsum, aes(x = new_chrom_end, xend = new_chrom_end, y = 0, yend = 10), lwd = 0.25, inherit.aes = FALSE, color = "grey75") +
  geom_segment(data = chrom_pos_cumsum, aes(x = new_chrom_end, xend = new_chrom_end, y = 0, yend = 7), lwd = 0.25, inherit.aes = FALSE, color = "grey75") +
  geom_segment(data = qtl_meta_use1_newpos, aes(x = new_qtl_pos - 1e6, xend = new_qtl_pos + 1e6, y = -0.5, yend = -0.5), lwd = 2, inherit.aes = FALSE) +
  geom_point(size = 0.1) +
  geom_point(data = gwas_mlmm_marker_info_ann_toplot, size = 0.8) + 
  # Known genes
  geom_text(data = known_genes_newpos1, aes(x = label_pos, y = 4, label = gene, fontface = "italic"), inherit.aes = FALSE, 
            check_overlap = TRUE, size = 2) + 
  geom_segment(data = known_genes_newpos1, aes(x = label_pos, xend = new_gene_pos, y = 3.5, yend = 2.5), lwd = 0.25, inherit.aes = FALSE) + 
  scale_color_manual(name = NULL, values = colors_use, guide = guide_legend(override.aes = list(size = 1.5))) +
  scale_shape_discrete(name = NULL) +
  # Remove annoying whitespace
  # scale_y_continuous(expand = c(0, 0)) + 
  scale_x_continuous(expand = c(0.005, 0), breaks = chrom_pos_cumsum$chrom_label_pos, labels = chrom_pos_cumsum$chrom) + 
  # ylab(expression(-log[10](italic(p)))) + 
  ylab(expression(-log[10](italic(q)))) + 
  ylim(c(-1, 8)) +
  # xlim(c(0, 5e9)) +
  facet_grid(trait ~ ., switch = "y", space = "free_x") + 
  theme_pnas() +
  theme(axis.ticks.x = element_blank(), axis.title.x = element_blank(), 
        legend.position = "bottom", legend.margin = margin(t = -5, unit = "pt"), legend.key.height = unit(0.75, "lines"),
        axis.line.x = element_blank(), strip.placement = "outside", panel.spacing.y = unit(0.5, "line"))



# Save
ggsave(filename = "gwas_manhattan_all_trait_QG_annotated.jpg", plot = g_all_gwas_annotated, path = fig_dir, 
       height = 8, width = 8.7, units = "cm", dpi = 1000)













### Test for pleiotropy

# Load the results
load(file.path(result_dir, "pheno_fw_gwas_pleiotropy_results.RData"))

## Tidy up the data
gwas_pleio_pheno_mean_fw_tidy <- gwas_pleio_pheno_mean_fw %>%
  filter(model == "QG") %>%
  gather(trait2, trait2_neg_log_p, b:log_delta) %>% 
  rename(trait1_neg_log_p = g) %>%
  filter(!is.na(trait2_neg_log_p),
         model == model_use) %>% 
  as_data_frame()

## The test for pleiotropy is a union-intersection test
## For each marker, there are two betas: beta_g and beta_stability,
## where beta_g is the marker effect for the genotype mean and beta_stability
## is the marker effect for phenotypic stability.
## 
## H0: beta_g == 0 or beta_stability == 0
## HA: beta_g != 0 and beta_stability != 0
## 
## We can use the maximum p-value as the test - if this p-value is less
## than a threshold alpha, the null hypothesis of no pleiotropy can be rejected.
## 

## Run a shortcut of the intersection-union test by grabbing the minimum -log10(p)
# Convert p-value to q-values
gwas_pleio_pheno_mean_fw_tidy_adj <- gwas_pleio_pheno_mean_fw_tidy %>% 
  mutate(min_neg_log_p = pmin(trait1_neg_log_p, trait2_neg_log_p),
         min_pvalue = 10^-min_neg_log_p) %>%
  group_by(trait, trait2) %>% 
  mutate(qvalue = qvalue(p = min_pvalue)$qvalue, 
         neg_log_q = -log10(qvalue),
         color = if_else(chrom %in% seq(1, 7, 2), "B", "G"),
         empirical_cutoff = quantile(min_neg_log_p, 1 - ((alpha / 10) / 2)),
         significant = min_neg_log_p > empirical_cutoff) %>%
  ungroup()


# Combine the minor allele frequency information
gwas_sig_pleio_marker_info <- gwas_pleio_pheno_mean_fw_tidy_adj %>% 
  group_by(trait, trait2) %>%
  top_n(n = 25, wt = min_neg_log_p) %>%
  ungroup() %>%
  left_join(., af1) %>%
  select(-af)

# Gather significant markers and create a GRange object
gwas_sig_pleio_grange <- gwas_sig_pleio_marker_info %>%
  # Add twice the window to complement the other overlapping procedures
  mutate(start = pos - (window), end = pos + (window)) %>%
  makeGRangesFromDataFrame(., keep.extra.columns = TRUE)

gwas_sig_pleio_grange_split <- gwas_sig_pleio_grange %>% 
  as.data.frame() %>%
  split(.$trait) %>%
  map(~split(., .$trait2) %>% map(~makeGRangesFromDataFrame(df = ., keep.extra.columns = TRUE)) %>% 
        GRangesList())




## Find the overlap with the MLMM markers
common_traits <- intersect(names(gwas_sig_pleio_grange_split), names(qtl_meta_grange))

gwas_sig_pleio_qtl_meta_overlap <- list(gwas_sig_pleio_grange_split[common_traits], qtl_meta_grange[common_traits]) %>%
  pmap(., ~lapply(as.list(.x), function(mar) mergeByOverlaps(query = mar, subject = .y))) %>%
  map(~do.call("rbind", .)) %>% do.call("rbind", .) %>%
  subset(x = ., , c(trait, trait2, marker, gene, reference)) %>%
  as_data_frame()


## Nest the marker overlaps and combine with the original data.frame for printing
gwas_sig_pleio_marker_info_ann <- gwas_sig_pleio_marker_info %>%
  left_join(., group_by(gwas_sig_pleio_qtl_meta_overlap, trait, trait2, marker) %>% 
              nest(marker:reference, .key = "annotation")) %>%
  mutate(qtl_hits = map_dbl(annotation, ~ifelse(is.null(.), 0, nrow(.))))

## Count the significant markers with overlaps
gwas_sig_pleio_marker_info_ann1 <- gwas_sig_pleio_marker_info_ann %>% 
  group_by(trait, trait2) %>% 
  summarize(prop_qtl_hits = mean(!map_lgl(annotation, is.null)), 
            sum_qtl_hits = sum(qtl_hits))

## Output a table
gwas_sig_pleio_marker_toprint <- gwas_sig_pleio_marker_info_ann %>%
  select(Trait = trait, Coef = trait2, Marker = marker, Chrom = chrom, Position = pos, 
         pvalue = min_pvalue, MAF = maf, QTLHits = qtl_hits, annotation) %>% 
  arrange(Trait, Coef, Chrom, Position) %>% 
  mutate(Coef = str_replace_all(Coef, coef_replace),
         Reference = map(annotation, "reference") %>% map(~ifelse(is.na(.), "T3", .)) %>% 
           map_chr(~paste0(unique(.), collapse = ", "))) %>%
  select(-annotation)

write_csv(x = gwas_sig_pleio_marker_toprint, path = file.path(fig_dir, "gwas_sig_plei_significant_associations.csv"))





# Manhattan plot
g_mod_man <- list(
  geom_point(),
  ylab(expression(-log[10](italic(p)))),
  xlab("Position (Mbp)"),
  scale_color_manual(values = color, guide = FALSE),
  theme_bw(),
  theme_manhattan(),
  theme(panel.border = element_blank())
)


# Labeller function
label_coef <- function(x) str_c("Genotype Mean\n", x)

## Manhattan plot for the pleitropic scan
## Iterate over each trait
## Add known QTL and color points based on those that overlap
for (tr in unique(gwas_pleio_pheno_mean_fw_tidy_adj$trait)) {
  
  gwas_pleio_toplot <- gwas_pleio_pheno_mean_fw_tidy_adj %>%
    filter(trait == tr) %>%
    left_join(., select(gwas_sig_pleio_marker_info_ann, trait, marker, trait2, min_neg_log_p, qtl_hits)) %>%
    mutate(qtl_hits = ifelse(is.na(qtl_hits), 0, qtl_hits),
           color = case_when(significant & qtl_hits == 0 ~ "Bl",
                             significant & qtl_hits > 0 ~ "Or",
                             TRUE ~ color),
           trait2 = str_replace_all(trait2, coef_replace))
                             
    
  g_gwas_plei <- gwas_pleio_toplot %>%
    ggplot(aes(x = pos / 1000000, y = min_neg_log_p, color = color)) + 
    g_mod_man +
    facet_grid(trait + trait2 ~ chrom, scales = "free", space = "free_x", switch = "x",
               labeller = labeller(trait2 = label_coef))
  
  
  save_file <- file.path(fig_dir, str_c("gwas_plei_manhattan_", tr, ".jpg"))
  ggsave(filename = save_file, plot = g_gwas_plei, width = 9, height = 6)
  
}


## Reformat the known genes
known_genes_newpos2 <- known_genes_newpos1 %>% 
  select(trait, trait:chrom, label_pos, new_gene_pos) %>%
  full_join(mutate(distinct(gwas_pleio_pheno_mean_fw_tidy_adj, trait, trait2), trait2 = str_replace_all(trait2, coef_replace))) %>%
  filter((gene == "Ppd-H2" & trait2 == "Linear Stability") | (gene == "Vrn-H1") | (gene == "Vrn-H3") )


## Output a pretty manhattan plot
g_gwas_pleio <- gwas_pleio_pheno_mean_fw_tidy_adj %>%
  left_join(snp_info_new_pos) %>%
  mutate(trait2 = str_replace_all(trait2, coef_replace),
         color = ifelse(min_neg_log_p >= empirical_cutoff, "Bl", color)) %>%
  ggplot(aes(x = new_pos, y = min_neg_log_p, color = color)) +
  # Add chromosome lines
  geom_segment(data = chrom_pos_cumsum, aes(x = new_chrom_end, xend = new_chrom_end, y = 0, yend = 5), lwd = 0.25, inherit.aes = FALSE, color = "grey75") +
  geom_point(size = 0.1) +
  geom_text(data = chrom_pos_cumsum, aes(x = chrom_label_pos, y = -1.5, label = chrom), size = 2, vjust = 0.70, inherit.aes = FALSE) +
  # Known genes
  geom_text(data = known_genes_newpos2, aes(x = label_pos, y = 4, label = gene, fontface = "italic"), inherit.aes = FALSE, 
            check_overlap = TRUE, size = 2) + 
  geom_segment(data = known_genes_newpos2, aes(x = label_pos, xend = new_gene_pos, y = 3.5, yend = 2.75), lwd = 0.25, inherit.aes = FALSE) +
  geom_segment(data = qtl_meta_use1_newpos, aes(x = new_qtl_pos - 1e6, xend = new_qtl_pos + 1e6, y = -0.5, yend = -0.5), lwd = 2, inherit.aes = FALSE) +
  scale_color_manual(name = NULL, values = color, guide = FALSE) +
  scale_x_continuous(expand = c(0.005, 0)) + 
  ylab(expression(-log[10](italic(p)))) + 
  ylim(c(-1.5, 5)) +
  facet_grid(trait + trait2 ~ ., switch = "y", space = "free_x", labeller = labeller(trait2 = function(x) str_c("Mean:", x))) + 
  theme_pnas() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), 
        legend.position = "bottom", legend.margin = margin(t = -5, unit = "pt"), legend.key.height = unit(0.75, "lines"),
        axis.line.x = element_blank(), strip.placement = "outside")
  
## Save
ggsave(filename = "gwas_pleio_manhattan_all_trait_QG_annotated.jpg", plot = g_gwas_pleio, path = fig_dir,
       height = 14, width = 11.4, units = "cm", dpi = 1000)
  











#### Resampling the Number of Environments

## Resample 20, 40, 60, or 80% of environments 250 times each, calculate the stability 
## estimates, then performed the mapping.

# Read in the results
load(file.path(result_dir, "pheno_fw_gwas_resample_results.RData"))


# Bind rows, unnest, and tidy
resample_gwas_sig_out1 <- resample_gwas_sig_out %>%
  bind_rows() %>%
  mutate(gwas_sig = map(results, "gwas_sig"),
         n_NA = map(results, "n_NA")) 


## Filter for the significant GWAS hits
resample_gwas_sig <- resample_gwas_sig_out1 %>%
  unnest(gwas_sig)

## Filter the original GWAS results for the stability QTL
gwas_sig_stability <- gwas_pheno_mean_fw_adj %>% 
  filter(coef != "g", qvalue <= alpha)

# Iterate over the resampling results and find the signficant loci at alpha
resample_gwas_sig_out2 <- resample_gwas_sig %>%
  ## Convert trait, p, iter, and coef to factors
  mutate_at(vars(trait, p, iter, coef), as.factor) %>%
  mutate(iter = factor(iter, levels = seq(250)))


## Since only the significant markers are part of this data.frame, in order
## to count the number of significant associations, we need to fill in unobserved
## combinations with 0
resample_gwas_sig_out_count <- resample_gwas_sig_out2 %>% 
  complete(trait, p, iter, coef) %>% 
  # First find the number of associations per experiment
  group_by(trait, p, coef, iter) %>%
  summarize(n_sig_SNP = sum(!is.na(marker))) %>%
  # Now take the mean and sd across iterations
  summarize_at(vars(n_sig_SNP), funs(mean = mean, sd = sd))


## Plot this
g_resample_n_SNPs <- resample_gwas_sig_out_count %>% 
  ggplot(aes(x = p, y = mean, col = trait, group = trait)) + 
  geom_point() + 
  geom_line() + 
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.25) + 
  ylab("Number of Significant SNPs") +
  xlab("Proportion of Environments Sampled") +
  facet_grid(~ coef) + 
  theme_bw()


## Find the intersection of the associations detected in the resampling experiment
## with the original associations.

# Subset the original gwas results
gwas_sig_fw_intersect <- gwas_sig_stability %>%
  select(trait, coef, marker)

# For each p, iter, trait, and coef, how often do we detect new SNPs or the same
# SNPs as in the full analysis?
resample_gwas_sig_fw_intersect <- resample_gwas_sig_out2 %>% 
  select(p, iter, marker, trait, coef) %>%
  group_by(p, iter, trait, coef) %>%
  do({
    df <- .
    
    # Trait and coef
    tr <- unique(df$trait)
    co <- unique(df$coef)
    
    # Filter the original data
    gwas_sig_fw_intersect1 <- gwas_sig_fw_intersect %>%
      filter(trait == tr, coef == co)
    
    # Intersect
    marker_intersect <- intersect(df$marker, gwas_sig_fw_intersect1$marker)
    marker_new <- setdiff(df$marker, gwas_sig_fw_intersect1$marker)
    
    # Return a data.frame of the intersected markers, the unique markers, and the counts of each
    data_frame(hit_SNPs = list(marker_intersect), 
               new_SNPs = list(marker_new))  })


## Again, we have to complete the observations, by filling in with 0s
# Then count the number of hit SNPs and new SNPS
resample_gwas_sig_fw_intersect_count <- resample_gwas_sig_fw_intersect %>% 
  ungroup() %>% 
  complete(trait, coef, p, iter) %>%
  mutate_at(vars(hit_SNPs, new_SNPs), ~map_int(., length)) %>%
  gather(SNP_type, count, hit_SNPs, new_SNPs)
  

## Summarize across iterations
resample_gwas_sig_fw_intersect_summary <- resample_gwas_sig_fw_intersect_count %>% 
  group_by(trait, coef, p, SNP_type) %>% 
  summarize_at(vars(count), funs(mean, sd)) %>%
  ungroup()


# Create a data.frame to plot
resample_gwas_sig_fw_intersect_toplot <- resample_gwas_sig_fw_intersect_summary %>% 
  mutate(coef = str_replace_all(coef, coef_replace), 
         SNP_type = str_replace_all(SNP_type, c("hit_SNPs" = "Recovered SNPs", "new_SNPs" = "New SNPs")))

# Plot
g_resample_intersect <- resample_gwas_sig_fw_intersect_toplot %>% 
  ggplot(aes(x = p, y = mean, fill = SNP_type, group = SNP_type)) + 
  geom_col(position = "dodge") + 
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.25, 
                position = position_dodge(0.9)) + 
  ylab("Number of SNPs") +
  xlab("Proportion of Environments Sampled") +
  scale_fill_discrete(name = "Association\nType") + 
  facet_grid(trait ~ coef) + 
  theme_bw() +
  theme(legend.position = "bottom")





## Count up the intersected and new SNPs and rank by number of times
## the SNP was detected

# First record the total number of iterations
n_iter <- nlevels(resample_gwas_sig_fw_intersect$iter)

resample_snp_detect_count <- resample_gwas_sig_fw_intersect %>% 
  ungroup() %>% 
  gather(SNP_type, marker, -p:-coef) %>%
  unnest() %>% 
  group_by(p, trait, coef, SNP_type, marker) %>% 
  summarize(times_detected = n(), prop_detected = times_detected / n_iter) %>%
  ungroup() %>%
  left_join(., select(S2TP_imputed_multi_genos_hmp, marker = rs, chrom, pos), by = "marker") %>%
  # Assign proportion bins 
  mutate(prop_group = cut(x = prop_detected, breaks = seq(0, 1, 0.2), include.lowest = TRUE, right = TRUE))

# Create a data.frame for plotting
resample_snp_detect_count_toplot <- resample_snp_detect_count %>%
  complete(p, trait, coef, SNP_type, fill = list(prop_detected = 0)) %>%
  mutate(coef = str_replace_all(coef, coef_replace), 
         SNP_type = str_replace_all(SNP_type, c("hit_SNPs" = "Recovered SNPs", "new_SNPs" = "New SNPs")))

# Are we more likely to discover new SNPs or those we discovered before?
g_resample_snp_detect_count <- resample_snp_detect_count_toplot %>%
  ggplot(aes(x = p, y = prop_detected, color = SNP_type)) +
  geom_boxplot(position = "dodge", width = 0.5) +
  scale_color_discrete(drop = FALSE) +
  ylim(c(0, 1))+
  facet_grid(trait ~ coef) + 
  ylab("Probability of Detecting a SNP") + 
  xlab("Proportion of Environments Sampled") +
  scale_color_discrete(name = "Association\nType") + 
  theme_bw() +
  theme(legend.position = "bottom")


## Combine the plots
# g_resample_combined <- g_resample_intersect + g_resample_snp_detect_count
g_resample_combined <- plot_grid(g_resample_intersect, g_resample_snp_detect_count, 
                                 labels = c("A", "B"))


save_file <- file.path(fig_dir, "resample_gwas_detection_results.jpg")
ggsave(filename = save_file, plot = g_resample_combined, width = 8, height = 5)





# Sort and display
resample_snp_detect_count %>%
  filter(times_detected == max(times_detected))


## Plot the detection results across chromosoems

# Create a modifier for the proportion bins
prop_group_replace <- expression(paste(0 <= x) <= 0.2, paste(0.2 < x) <= 0.4, paste(0.4 < x) <= 0.6, paste(0.6 < x) <= 0.8, paste(0.8 < x) <= 1)


# Common plot modifier
g_mod <- list(geom_point(),
              geom_segment(data = barley_lengths, aes(x = 0, y = 0, xend = 0, yend = length / 1000000), 
                           inherit.aes = FALSE, size = 2, lineend = "round", col = "grey70") ,
              ylab("Position (Mbp)") ,
              xlab("Chromosome"),
              xlim(c(-2, 2)) , 
              facet_grid(~ chrom, switch = "x") , 
              scale_shape_manual(values = c("Recovered SNPs" = 1, "New SNPs" = 0), name = "Association\nType"), 
              scale_size_manual(name = "Probability of\nDetecting SNP", labels = prop_group_replace,
                                values = unique(as.numeric(resample_snp_detect_count$prop_group))),                  
              scale_color_discrete(name = "Stability\nType"),
              scale_y_reverse() , 
              theme_bw() , 
              theme(panel.grid = element_blank(), 
                    panel.border = element_blank(), 
                    # panel.spacing.x = unit(1, units = "cm"),
                    axis.text.x = element_blank(), 
                    axis.title.x = element_blank(), 
                    axis.ticks.x = element_blank(),
                    strip.background = element_blank(),
                    legend.justification = c(0.5, 1)) # Set the legend justification in the plot area
)


# Create a data.frame to plot
resample_snp_detect_count_toplot <- resample_snp_detect_count %>%
  mutate(p = parse_number(p),
         # Adjust the x axis position to "jitter" the points
         x = if_else(SNP_type == "hit_SNPs", -p * 2, p * 2),
         coef = str_replace_all(coef, coef_replace),
         SNP_type = str_replace_all(SNP_type, c("hit_SNPs" = "Recovered SNPs", "new_SNPs" = "New SNPs")))


## Create a legend for the positioning of different environment proportion in relation to
## the chromosome
g_resample_legend <- distinct(resample_snp_detect_count_toplot, p, x, SNP_type) %>%
  ggplot(aes(x = x, y = 3, shape = SNP_type)) +
  geom_point(size = 2) +
  geom_text(aes(x = x, y = 2, label = p), angle = 45, size = 2.5) + 
  geom_segment(aes(x = 0, y = 1, xend = 0, yend = 5), 
               inherit.aes = FALSE, size = 2, lineend = "round", col = "grey70") + 
  scale_shape_manual(values = c("Recovered SNPs" = 1, "New SNPs" = 0), guide = FALSE) +
  xlab("Proportion of\nEnvironments Sampled") +
  labs(title = "Proportion of\nEnvironments Sampled") + 
  theme_bw() +
  theme(panel.grid = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        title = element_text(size = 8),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.margin = unit(c(0,0,0,0), "cm"))

# Save
save_file <- file.path(fig_dir, "resample_legend.jpg")
ggsave(filename = save_file, plot = g_resample_legend, width = 1.5, height = 1.5, dpi = 1000)




# Create a 'heatmap' for each trait
g_resample_hd <- resample_snp_detect_count_toplot %>% 
  filter(trait == "HeadingDate") %>% 
  ggplot(aes(y = pos / 1000000, x = x, size = prop_group, col = coef, shape = SNP_type)) + 
  g_mod +
  labs(title = "Heading Date")


g_resample_ph <- resample_snp_detect_count_toplot %>% 
  filter(trait == "PlantHeight") %>% 
  ggplot(aes(y = pos / 1000000, x = x, size = prop_group, col = coef, shape = SNP_type)) + 
  g_mod +
  labs(title = "Plant Height")

g_resample_gy <- resample_snp_detect_count_toplot %>% 
  filter(trait == "GrainYield") %>% 
  ggplot(aes(y = pos / 1000000, x = x, size = prop_group, col = coef, shape = SNP_type)) + 
  g_mod +
  labs(title = "Grain Yield")

## Plot together
save_file <- file.path(fig_dir, "resample_count_hd.jpg")
ggsave(filename = save_file, plot = g_resample_hd, width = 10, height = 7, dpi = 1000)

save_file <- file.path(fig_dir, "resample_count_ph.jpg")
ggsave(filename = save_file, plot = g_resample_ph, width = 10, height = 7, dpi = 1000)

save_file <- file.path(fig_dir, "resample_count_gy.jpg")
ggsave(filename = save_file, plot = g_resample_gy, width = 10, height = 7, dpi = 1000)




