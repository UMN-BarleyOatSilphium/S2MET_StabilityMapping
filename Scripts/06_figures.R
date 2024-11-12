## S2MET Stability Mapping
## 
## figures.R
## 
## Script for generating figures
## 



# Project and other directories
source(here::here("startup.R"))

library(patchwork)
library(cowplot)

# Set a window outside of each significant SNP to search for genes
sig_snp_window <- 4.6 * 1e6

# New colors for different parameters
parameter_colors_use <- set_names(neyhart_palette("umn1")[3:5], names(coef_replace))

# Read in known genes
barley_known_genes <- read_csv(file = file.path(data_dir, "known_genes_qtl.csv")) %>%
  mutate(chrom1 = paste0("chr", chrom))

# Colors for chromosomes
# chrom_colors <- neyhart_palette("umn2")[1:7]
chrom_colors <- c("#d95218", "#e5a000", "#7f2e8d", "#5ea500", "#0095d5", "#a31431", "#0c55a9")
names(chrom_colors) <- paste0("chr", seq_along(chrom_colors))




# Ad hoc functions --------------------------------------------------------

# Function to display an LD plot
ld_plot <- function (geno, map, by.chrom = FALSE, max.pair = 10000, dof = 8, max.loci = NULL, position = "bp") {
  require(scam)
  chroms <- levels(as.factor(map$chrom))
  n.chrom <- length(chroms)
  result <- NULL
  for (i in 1:n.chrom) {
    ix <- which(map$chrom == chroms[i])
    m <- length(ix)
    if (!is.null(max.loci)) {
      if (m > max.loci) {
        ix <- sample(ix, max.loci)
        m <- max.loci
      }
    }
    tmp <- expand.grid(col = 1:m, row = 1:m)
    tmp <- tmp[tmp$row >= tmp$col, ]
    r2 <- cor(geno[, ix])^2
    r2.vec <- as.vector(r2[cbind(tmp$row, tmp$col)])
    d <- as.matrix(dist(matrix(map$pos[ix], ncol = 1)))
    d.vec <- as.vector(d[cbind(tmp$row, tmp$col)])
    if (position == "bp") 
      d.vec <- d.vec/1e+06
    result <- rbind(result, data.frame(chrom = chroms[i], d = d.vec, r2 = r2.vec))
  }
  if (by.chrom) {
    spline.data <- NULL
    for (i in 1:n.chrom) {
      result_i <- result[result$chrom == chroms[i],]
      max.pair <- min(max.pair, nrow(result_i))
      scam.ans <- scam(formula = r2 ~ s(d, bs = c("mdcx"), k = dof), 
                       data = result_i[sample(1:nrow(result_i), max.pair), ])
      dmax <- max(result_i$d)
      predans <- predict.scam(scam.ans, newdata = data.frame(d = seq(0,  dmax, length.out = 500)))
      spline.data <- rbind(spline.data, data.frame(chrom = chroms[i], d = seq(0, dmax, length.out = 500), r2 = predans))
    }
    p <- ggplot(data = spline.data, aes(x = d, y = r2, color = chrom)) + 
      ylab(expression(r^2)) + 
      theme_bw() + 
      geom_line()
  } else {
    max.pair <- min(max.pair, nrow(result))
    scam.ans <- scam(formula = r2 ~ s(d, bs = c("mdcx"), k = dof), 
                     data = result[sample(1:nrow(result), max.pair), ])
    dmax <- max(result$d)
    predans <- predict.scam(scam.ans, newdata = data.frame(d = seq(0,  dmax, length.out = 500)))
    spline.data <- data.frame(d = seq(0, dmax, length.out = 500), r2 = predans)
    p <- ggplot(data = spline.data, aes(x = d, y = r2)) + 
      ylab(expression(r^2)) + 
      theme_bw() + geom_line()
  }
  if (position == "bp") {
    p <- p + xlab("Distance (Mb)")
  }
  else {
    p <- p + xlab("Distance (cM)")
  }
  return(p)
}






# Load results ------------------------------------------------------------

load(file.path(data_dir, "project_data.RData"))
load(file.path(result_dir, "stability_estimates.RData"))
load(file.path(result_dir, "stability_gwas_out.RData"))
load(file.path(result_dir, "univariate_gwas_analysis.RData"))
load(file.path(result_dir, "genomewide_prediction_results.RData"))




# Figure 1: Reaction norms and phenotypic correlations ----------

# Plot reaction norms
g_reaction_list <- stability_data_use %>%
  filter(line_name %in% tp, trait %in% traits) %>%
  group_by(trait) %>%
  do(plot = {
    df <- .
    
    # Subset the stability estimates
    stab_data <- genotype_stability %>%
      filter(trait == unique(df$trait)) %>%
      select(-ranova) %>%
      unnest(regression) %>%
      filter(line_name %in% tp)
    
    ggplot(data = df, aes(x = h, y = value, group = line_name)) +
      geom_point() +
      # geom_smooth(method = "lm", se = FALSE, formula = y ~ x) +
      geom_abline(data = stab_data, aes(intercept = g, slope = b, color = b), alpha = 0.75) +
      scale_color_gradient2(midpoint = 1, name = expression(italic(b[i])), guide = guide_colorbar(title.position = "left")) +
      scale_x_continuous(name = expression('Environmental effect'~(italic(t[j]))), breaks = pretty) +
      scale_y_continuous(name = "Phenotypic value", breaks = pretty) +
      theme_presentation() +
      theme(legend.position = "inside", legend.position.inside = c(0.05, 0.95), legend.justification = c(0, 1))
    
   }) %>% ungroup()


# Combine
g_reaction_combined <- stability_data_use %>%
  filter(line_name %in% tp, trait %in% traits) %>%
  left_join(., unnest(select(genotype_stability, -ranova), regression)) %>%
  mutate(trait = paste0("'", str_add_space(trait), " ('*", str_replace_all(trait_units[trait], " ", "~"), "*')'")) %>%
  ggplot(aes(x = h, y = value, group = line_name)) +
  geom_point(size = 0.5) +
  # geom_smooth(method = "lm", se = FALSE, formula = y ~ x) +
  geom_abline(aes(intercept = g, slope = b, color = b), alpha = 0.75, linewidth = 0.5) +
  scale_color_gradient2(midpoint = 1, name = expression(italic(b[i])), guide = guide_colorbar(title.position = "left")) +
  scale_x_continuous(name = expression('Environmental effect'~(italic(t[j]))), breaks = pretty) +
  scale_y_continuous(name = "Phenotypic value", breaks = pretty) +
  facet_wrap(~ trait, ncol = 1, scales = "free", strip.position = "left", labeller = labeller(trait = label_parsed)) +
  theme_presentation(base_size = 10) +
  theme(legend.position = "none", strip.placement = "outside")

# Combine mean-stability plot
genotype_stability1 <- genotype_stability %>%
  select(-ranova) %>%
  unnest(regression) %>%
  filter(trait %in% traits, line_name %in% tp) %>%
  mutate(sd_d = sqrt(var_d))

genotype_stability_ann <- genotype_stability1 %>%
  group_by(trait) %>%
  summarize(xAnnot = max(g) - (0.15 * diff(range(g))), yAnnot = max(b) - (0.95 * diff(range(b))),
            corr = cor(g, b), .groups = "drop") %>%
  mutate(annot = paste0("italic(r)==", format_numbers(corr)))

g_mean_linStability <- genotype_stability1 %>%
  ggplot(data = , mapping = aes(x = g, y = b)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
  geom_text(data = genotype_stability_ann, aes(x = xAnnot, y = yAnnot, label = annot), parse = TRUE, size = 3) +
  scale_x_continuous(name = expression('Genotype mean'~(italic(g[i]))), breaks = pretty) +
  scale_y_continuous(name = expression('Slope'~(italic(b[i]))), breaks = pretty) +
  facet_wrap(~ trait, ncol = 1, scales = "free", strip.position = "left", labeller = labeller(trait = label_parsed)) +
  theme_presentation(base_size = 10)
  

genotype_stability_ann <- genotype_stability1 %>%
  group_by(trait) %>%
  summarize(xAnnot = max(g) - (0.15 * diff(range(g))), yAnnot = max(sd_d) - (0.95 * diff(range(sd_d))),
            corr = cor(g, sd_d), .groups = "drop") %>%
  mutate(annot = paste0("italic(r)==", format_numbers(corr)))

g_mean_nonLinStability <- genotype_stability1 %>%
  ggplot(data = , mapping = aes(x = g, y = sd_d)) +
  geom_point(size = 0.5) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
  geom_text(data = genotype_stability_ann, aes(x = xAnnot, y = yAnnot, label = annot), parse = TRUE, size = 3) +
  scale_x_continuous(name = expression('Genotype mean'~(italic(g[i]))), breaks = pretty) +
  scale_y_continuous(name = expression('Deviation '~(italic(sigma[delta(i)]))), breaks = pretty) +
  facet_wrap(~ trait, ncol = 1, scales = "free", strip.position = "left", labeller = labeller(trait = label_parsed)) +
  theme_presentation(base_size = 10)



# Combine the plots
g_combined <- wrap_plots(g_reaction_combined, 
                         g_mean_linStability + theme(strip.text = element_blank()), 
                         g_mean_nonLinStability + theme(strip.text = element_blank())) + 
  plot_annotation(tag_levels = "A")

# Save
ggsave(filename = "figure1_reaction_stability_correlations.png", plot = g_combined, path = fig_dir,
       width = 7, height = 8, dpi = 500)







# Figure 2: combined GWAS manhattan plot ----------------------------------

# Set the significance threshold by calculate the effective number of markers
# Fist calculate the squared correlations between all markers in each chromosome
me_list <- snp_info %>%
  split(.$chrom) %>%
  map_dbl(~{
    mat1 <- geno_mat_mapping[,.x$marker]
    r2 <- cor(mat1)^2
    GWASpoly:::Keff(r2 = r2, alpha = 0.05)
  })
me <- sum(me_list)
thresh <- -log10(0.05 / me)





# Add the position of the preceding chromosome to the position of the next
chrom_pos_cumsum <- snp_info %>% 
  group_by(chrom) %>% 
  summarize(chrom_pos = max(pos)) %>% 
  mutate(max_pos = c(0, head(chrom_pos, -1)), 
         next_chrom_start_pos = cumsum(max_pos),
         new_chrom_end = chrom_pos + next_chrom_start_pos,
         chrom_label_pos =( next_chrom_start_pos + new_chrom_end) / 2)

# Create a list of SNPs, old position, and new position
snp_info_plot <- snp_info %>%
  left_join(., chrom_pos_cumsum) %>%
  mutate(new_pos = pos + next_chrom_start_pos) %>%
  select(marker, chrom, pos, new_pos)


## Supplemental plots - manhattan plots of each trait / parameter

gwas_toplot <- univariate_gwas_farm_bestModel_qq %>%
  filter(trait %in% traits) %>%
  select(trait, parameter, nPC = bestModel_qq_nPC) %>%
  left_join(., univariate_gwas_farm) %>%
  select(trait, parameter, nPC, scores, fdr_p_thresh) %>%
  mutate(scores = map(scores, ~merge(.x, snp_info_plot) %>% arrange(marker)))

# DF of thresholds
gwas_thresh <- gwas_toplot %>%
  distinct(trait, parameter, fdr_p_thresh) %>%
  mutate(thresh_score = -log10(fdr_p_thresh))


# Plot
gwas_plot_list <- gwas_toplot %>%
  group_by(trait) %>%
  nest() %>%
  ungroup() %>%
  mutate(plot = list(NULL))

for (i in seq_len(nrow(gwas_plot_list))) {
  
  df <- gwas_plot_list$data[[i]]
  
  scores <- df %>%
    unnest(scores) %>%
    mutate(score = -log10(P.value),
           fdr_p_thresh_use = ifelse(is.na(fdr_p_thresh), -Inf, fdr_p_thresh))
  
  chrom_pos_cumsum1 <- chrom_pos_cumsum %>%
    mutate(yend = ceiling(1.1 * max(scores$score)))
    
  gplot <- ggplot(data = scores, aes(x = new_pos, y = score)) + 
    # Add chromosome breakpoints
    geom_segment(data = chrom_pos_cumsum1, aes(x = new_chrom_end, xend = new_chrom_end, y = 0, yend = yend), 
                 inherit.aes = FALSE, color = "grey85") +
    # Add points, large points = significant markers
    geom_point(size = ifelse(scores$P.value <= scores$fdr_p_thresh_use, 2, 1)) +
    # Add chromosome names
    geom_text(data = chrom_pos_cumsum1, aes(x = chrom_label_pos, y = -1, label = chrom), size = 2, vjust = 0, inherit.aes = FALSE) +
    facet_grid(parameter + nPC ~ ., space = "free_x", scales = "free_x", switch = "y", 
               labeller = labeller(parameter = as_mapper(~coef_replace[.]), nPC = label_both)) +
    scale_x_continuous(name = "Chromosome", expand = expansion()) +
    scale_y_continuous(name = expression(-log[10](italic(p))), breaks = pretty, expand = expansion(mult = c(0, .1))) +
    labs(title = str_add_space(gwas_plot_list$trait[i])) +
    theme_genetics() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank(), 
          axis.line.x = element_blank(), legend.position = "bottom", legend.text = element_text(size = 6))
  
  gwas_plot_list$plot[[i]] <- gplot
  
}


# Save plots
for (i in seq_len(nrow(gwas_plot_list))) {
  trt <- gwas_plot_list$trait[i]
  plt <- gwas_plot_list$plot[[i]]
  
  filename <- paste0("univariate_gwas_manhattan_", trt, ".png")
  ggsave(filename = filename, plot = plt, path = fig_dir, width = 10, height = 7)
  
}



# Create combined manhattan plots with gene/QTL annotations


# Create a data frame of T3 QTL for these traits
t3_qtl_positions <- t3_qtl_granges %>%
  subset(t3_trait %in% c("grain protein", "grain yield", "heading date", "heading date (Julian)", "plant height", "test weight")) %>%
  as.data.frame(x = ., row.names = NULL, stringsAsFactors = FALSE) %>%
  select(-attributes) %>%
  mutate(trait = case_when(
    t3_trait %in% c("diastatic power", "grain protein", "alpha amylase") ~ "GrainProtein",
    t3_trait %in% c("grain yield") ~ "GrainYield",
    t3_trait %in% c("heading date", "heading date (Julian)") ~ "HeadingDate",
    t3_trait %in% c("plant height") ~ "PlantHeight",
    t3_trait %in% c("plump grain", "breeders plump grain") ~ "PlumpGrain",
    t3_trait %in% c("test weight") ~ "TestWeight"
  )) %>%
  # Distinctify
  distinct(seqnames, start, t3_trait, trait, phenotyping_source)


qtl_meta_use1_newpos <- t3_qtl_positions %>%
  mutate(chrom = parse_number(as.character(seqnames))) %>%
  left_join(., chrom_pos_cumsum) %>%
  mutate(new_qtl_pos = start + next_chrom_start_pos,
         chrom1 = as.factor(paste0("chr", chrom)))

gwas_toplot1 <- gwas_toplot %>%
  mutate(fdr_p_thresh = ifelse(is.na(fdr_p_thresh), -Inf, fdr_p_thresh),
         scores = map2(scores, fdr_p_thresh, ~mutate(.x, significant = ifelse(P.value <= .y, "sig", "nonsig"), score = -log10(P.value)))) %>%
  select(-fdr_p_thresh) %>%
  unnest(scores) %>%
  mutate(parameter = factor(parameter, levels = (names(coef_replace))))

names(univariate_gwas_overlap_perTrialT3QTL)[grepl("t3_qtl_granges", names(univariate_gwas_overlap_perTrialT3QTL))] <- "t3_qtl_granges"
names(univariate_gwas_overlap_experimentSetT3QTL)[grepl("t3_qtl_granges", names(univariate_gwas_overlap_experimentSetT3QTL))] <- "t3_qtl_granges"



## Add color for T3 QTL that overlap 
t3_qtl_perTrial_overlap_positions <- univariate_gwas_overlap_perTrialT3QTL$t3_qtl_granges %>%
  as.data.frame(x = ., row.names = NULL, stringsAsFactors = FALSE) %>%
  select(-attributes) %>%
  # Add marker, trait, and parameter
  mutate(trait = univariate_gwas_overlap_perTrialT3QTL$trait, parameter = univariate_gwas_overlap_perTrialT3QTL$parameter,
         marker = univariate_gwas_overlap_perTrialT3QTL$marker, marker_pos = start(ranges(univariate_gwas_overlap_perTrialT3QTL$univariate_gwas_sigmar_grange))) %>%
  # Subset traits
  inner_join(., dplyr::rename(distinct(qtl_meta_use1_newpos, seqnames, chrom, start, t3_trait, trait, new_qtl_pos), trait1 = trait)) %>%
  filter(trait == trait1) %>%
  mutate(parameter = factor(parameter, levels = names(coef_replace)),
         chrom1 = as.factor(paste0("chr", chrom))) %>%
  subset(parameter != "sd_d")

## Add color for T3 QTL that overlap 
t3_qtl_experiment_overlap_positions <- univariate_gwas_overlap_experimentSetT3QTL$t3_qtl_granges %>%
  as.data.frame(x = ., row.names = NULL, stringsAsFactors = FALSE) %>%
  select(-attributes) %>%
  # Add marker, trait, and parameter
  mutate(trait = univariate_gwas_overlap_experimentSetT3QTL$trait, parameter = univariate_gwas_overlap_experimentSetT3QTL$parameter,
         marker = univariate_gwas_overlap_experimentSetT3QTL$marker, marker_pos = start(ranges(univariate_gwas_overlap_experimentSetT3QTL$univariate_gwas_sigmar_grange))) %>%
  # Subset traits
  inner_join(., dplyr::rename(distinct(qtl_meta_use1_newpos, seqnames, chrom, start, t3_trait, trait, new_qtl_pos), trait1 = trait)) %>%
  filter(trait == trait1)

# Add named gene annotations
gene_overlap_positions <- univariate_gwas_sigmar_grange_geneOverlap %>%
  subset(!is.na(Name), barley_txdb2_genes_grange, drop = TRUE) %>%
  as.data.frame(x = ., row.names = NULL, stringsAsFactors = FALSE) %>%
  mutate(pos = (start + end) / 2,
         chrom = parse_number(as.character(seqnames)),
         chrom1 = as.factor(paste0("chr", chrom))) %>%
  left_join(., chrom_pos_cumsum) %>%
  mutate(new_pos = pos + next_chrom_start_pos)


# Remove non-linear stability
gwas_toplot2 <- gwas_toplot1 %>%
  # subset(parameter != "sd_d") %>%
  mutate(chrom1 = as.factor(paste0("chr", chrom)))

# Identify genes that are within the window for significant SNPs
genes_to_plot <- gwas_toplot2 %>%
  subset(significant == "sig") %>%
  mutate(.id = seq_len(nrow(.))) %>%
  split(.$.id) %>%
  map_df(~{
    genes1 <- subset(barley_known_genes, chrom == .x$chrom & start >= .x$pos - sig_snp_window & start <= .x$pos + sig_snp_window)
    genes1 <- merge(.x[,c("trait", "parameter", "chrom1", "score")], distinct(genes1, gene, chrom, chrom1, start)) %>%
      dplyr::rename(pos = start)
    genes2 <- subset(gene_overlap_positions, chrom == .x$chrom & pos >= .x$pos - sig_snp_window & pos <= .x$pos + sig_snp_window)
    genes2 <- genes2 %>%
      distinct() %>%
      select(chrom, chrom1, pos, gene = Name)
    genes2 <- merge(.x[,c("trait", "parameter", "chrom1", "score")], genes2)
    if (nrow(genes1) > 0) {
      top_n(genes1, n = 1, wt = -abs(pos - .x$pos))
    } else {
      genes_print <- bind_rows(genes1, genes2) %>%
        top_n(n = 1, wt = -abs(pos - .x$pos))
      genes_print
    }
  })
genes_to_plot1 <- genes_to_plot %>%
  subset(nchar(gene) < 10) %>%
  mutate(y1 = ceiling(score),
         y2 = max(y1),
         y_nudge = y2 - y1,
         x_nudge = ifelse(pos > 300e6, -100e6, 100e6)) %>%
  group_by(chrom1, trait, parameter, gene) %>%
  top_n(., n = 1, wt = y1) %>%
  ungroup() %>%
  select(-score)


## All traits and genes

g_gwas_unified <- gwas_toplot2 %>%
  ggplot(aes(x = pos, y = score, color = chrom1)) +
  # geom_point() +
  geom_segment(aes(y = 0, yend = score)) +
  geom_hline(data = gwas_thresh, aes(yintercept = thresh_score), lwd = 0.25) +
  # Add T3 QTL annotation
  geom_segment(data = qtl_meta_use1_newpos, aes(x = start - 1e5, xend = start + 1e5, y = -0.5, yend = -0.5, lty = "T3 QTL"),
               linewidth = 2, inherit.aes = FALSE) +
  geom_segment(data = t3_qtl_perTrial_overlap_positions, aes(x = start - 5e6, xend = start + 5e6, y = -0.5, yend = -0.5),
               linewidth = 2) +
  ggrepel::geom_text_repel(data = genes_to_plot1, aes(x = pos, y = y1, label = gene), size = 2, angle = 90, vjust = 0, hjust = 0, nudge_y = genes_to_plot1$y_nudge,
                           nudge_x = genes_to_plot1$x_nudge, box.padding = unit(0.2, "lines"), direction = "x", min.segment.length = 0.1, color = "black") +
  scale_x_continuous(name = "Chromosome", expand = expansion(), breaks = NULL, labels = NULL) +
  scale_y_continuous(name = expression(-log[10](italic(p))), breaks = pretty, expand = expansion(mult = c(0, .1))) +
  scale_color_manual(values = chrom_colors) +
  scale_linetype_discrete(name = NULL, guide = guide_legend(keywidth = unit(0.2, "line"))) +
  facet_grid(trait + parameter ~ chrom, scales = "free_x", space = "free_x", switch = "both", 
             labeller = labeller(trait = function(x) abbreviate(x, 2), parameter = coef_replace1)) +
  theme_classic() +
  theme(strip.background = element_blank(), panel.spacing.x = unit(0, "lines"), legend.position = "none",
        strip.placement = "outside", strip.clip = "off")


# Save
ggsave(filename = "figureS2_combined_gwas_results.png", plot = g_gwas_unified, path = fig_dir, 
       height = 12, width = 6, dpi = 1000)

# Subset select traits
traits_select <- c("GrainProtein", "GrainYield", "HeadingDate")
genes_to_plot2 <- subset(genes_to_plot1, trait %in% traits_select)
gwas_thresh2 <- subset(gwas_thresh, trait %in% traits_select & parameter != "sd_d")

g_gwas_unified_select <- gwas_toplot2 %>%
  subset(trait %in% traits_select & parameter != "sd_d") %>%
  ggplot(aes(x = pos, y = score, color = chrom1)) +
  # geom_point() +
  geom_segment(aes(y = 0, yend = score)) +
  geom_hline(data = gwas_thresh2, aes(yintercept = thresh_score), lwd = 0.25) +
  # Add T3 QTL annotation
  geom_segment(data = subset(qtl_meta_use1_newpos, trait %in% traits_select), aes(x = start - 1e5, xend = start + 1e5, y = -0.5, yend = -0.5, lty = "T3 QTL"),
               linewidth = 2, inherit.aes = FALSE) +
  geom_segment(data = subset(t3_qtl_perTrial_overlap_positions, trait %in% traits_select), aes(x = start - 5e6, xend = start + 5e6, y = -0.5, yend = -0.5),
               linewidth = 2) +
  ggrepel::geom_text_repel(data = genes_to_plot2, aes(x = pos, y = y1, label = gene), size = 2, angle = 90, vjust = 0, hjust = 0, nudge_y = genes_to_plot2$y_nudge,
                           nudge_x = genes_to_plot2$x_nudge, box.padding = unit(0.2, "lines"), direction = "x", min.segment.length = 0.1, color = "black") +  
  scale_x_continuous(name = "Chromosome", expand = expansion(), breaks = NULL, labels = NULL) +
  scale_y_continuous(name = expression(-log[10](italic(p))), breaks = pretty, expand = expansion(mult = c(0, .1))) +
  scale_color_manual(values = chrom_colors) +
  scale_linetype_discrete(name = NULL, guide = guide_legend(keywidth = unit(0.2, "line"))) +
  facet_grid(trait + parameter ~ chrom, scales = "free_x", space = "free_x", switch = "both", 
             labeller = labeller(trait = function(x) abbreviate(x, 2), parameter = coef_replace1, .multi_line = FALSE)) +
  theme_classic() +
  theme(strip.background = element_blank(), panel.spacing.x = unit(0, "lines"), legend.position = "none",
        strip.placement = "outside", strip.clip = "off")

# Save
ggsave(filename = "figure2_combined_gwas_results.png", plot = g_gwas_unified_select, path = fig_dir, 
       height = 6, width = 5, dpi = 1000)




## Calculate the distance between each significant association and the nearest 
## T3 QTL
gwas_t3_qtl_closest_marker <- t3_qtl_overlap_positions %>%
  select(chrom, start, t3_trait, trait, trait1, parameter, marker, marker_pos) %>%
  distinct() %>%
  mutate(marker_qtl_distance = abs(marker_pos - start)) %>%
  group_by(trait, parameter, marker) %>%
  top_n(x = ., n = 1, wt = -marker_qtl_distance) %>%
  ungroup()

# Calculate mean, min, max
gwas_t3_qtl_closest_marker %>%
  group_by(trait, parameter) %>%
  summarize_at(vars(marker_qtl_distance), list(mean = mean, min = min, max = max))


# trait        parameter     mean    min     max
# 1 GrainProtein b          574165.  52091 1998371
# 2 GrainProtein g          564718.  54482 3044882
# 3 GrainYield   b          286499   37095  989908
# 4 GrainYield   g          517494.    336 2374806
# 5 HeadingDate  b          438924.    565 2565823
# 6 HeadingDate  g          339734.    339  899792
# 7 PlantHeight  b          363538.    267  896281
# 8 PlantHeight  g         1171860   24906 2146308
# 9 PlumpGrain   b          506848.  20773 1856067
# 10 PlumpGrain   g          515802.   9011 2165569
# 11 PlumpGrain   sd_d       384870. 281162  553082
# 12 TestWeight   b          221925  221925  221925
# 13 TestWeight   g          631986     445 1275269



# Do the same thing for nearby genes
gwas_gene_overlap <- univariate_gwas_sigmar_grange_geneOverlap$barley_txdb2_genes_grange %>%
  as.data.frame(x = ., row.names = NULL, stringsAsFactors = FALSE) %>%
  # Add marker, trait, and parameter
  mutate(trait = univariate_gwas_sigmar_grange_geneOverlap$trait, parameter = univariate_gwas_sigmar_grange_geneOverlap$parameter,
         marker = univariate_gwas_sigmar_grange_geneOverlap$marker, marker_pos = start(ranges(univariate_gwas_sigmar_grange_geneOverlap$univariate_gwas_sigmar_grange)))


gwas_gene_overlap_closest_marker <- gwas_gene_overlap %>%
  select(start, end, gene_id, Name, trait, parameter, marker, marker_pos) %>%
  distinct() %>%
  mutate(marker_gene_distance = pmin(abs(marker_pos - start), abs(marker_pos - end))) %>%
  group_by(trait, parameter, marker) %>%
  top_n(x = ., n = 1, wt = -marker_gene_distance) %>%
  ungroup()

# Calculate mean, min, max
gwas_gene_overlap_closest_marker %>%
  group_by(trait, parameter) %>%
  summarize_at(vars(marker_gene_distance), list(mean = mean, min = min, max = max))


# trait        parameter   mean   min    max
# 1 GrainProtein b         15690.   136  71397
# 2 GrainProtein g         17536.   330  71397
# 3 GrainYield   b         10581.   180  27847
# 4 GrainYield   g         17606     29  86544
# 5 HeadingDate  b          4785.   201  33555
# 6 HeadingDate  g         18740    161 151249
# 7 PlantHeight  b         13474.   167  85530
# 8 PlantHeight  g          5621.   720  10344
# 9 PlumpGrain   b          5906.   532  12140
# 10 PlumpGrain   g          4833.    14  25529
# 11 PlumpGrain   sd_d       7348.    64  24336
# 12 TestWeight   b           281    281    281
# 13 TestWeight   g          2761     78  10064


gwas_gene_overlap %>%
  select(start, end, gene_id, Name, trait, parameter, marker, marker_pos) %>%
  distinct() %>%
  mutate(marker_gene_distance = pmin(abs(marker_pos - start), abs(marker_pos - end))) %>% 
  filter(!is.na(Name)) %>% 
  distinct(trait, parameter, marker, marker_gene_distance, Name, gene_id)










# Figure 3: Results of resampling -----------------------------------------


# Plot the mean-absolute-deviation of the sampled estimates from the full-data estimates

stability_sample_mae <- stability_sample_deviations %>%
  group_by(trait, pEnv, nEnv, rep) %>%
  summarize_at(vars(contains("dev")), list(mae = mean), na.rm = TRUE) %>%
  ungroup()

# Calculate the mean and 95% interval for each parameter
stability_sample_mae_toplot <- stability_sample_mae %>%
  dplyr::rename(sd_d_sample_dev_mae = d_sample_dev_mae) %>%
  gather(parameter, sample_mae, contains("mae")) %>%
  mutate(parameter = str_remove(parameter, "_sample_dev_mae")) %>%
  group_by(trait, parameter, pEnv) %>%
  summarize_at(vars(sample_mae), list(mean = mean, lower = ~quantile(., 0.05/2), upper = ~quantile(., 0.975))) %>%
  ungroup() %>%
  mutate(parameter = factor(parameter, levels = names(coef_replace1))) %>%
  subset(trait %in% traits)


# Split by trait and plot
stability_mae_plotlist <- stability_sample_mae_toplot %>%
  filter(trait %in% traits) %>%
  split(list(.$parameter, .$trait)) %>%
  map(~{
    # Plot
    g <- ggplot(.x, aes(x = pEnv, y = mean, color = parameter)) +
      geom_ribbon(aes(ymin = lower, ymax = upper, fill = parameter), alpha = 0.2) +
      geom_line() +
      # geom_point() +
      scale_x_continuous(name = "Proportion of total environments", limits = c(0.1, 0.9), breaks = seq(0.2, 0.8, by = 0.2)) +
      scale_y_continuous(name = 'Mean absolute deviation', breaks = pretty, limits = c(0, NA)) + 
      scale_color_manual(values = parameter_colors_use, labels = coef_replace1, name = NULL, guide = "none") +
      scale_fill_manual(values = parameter_colors_use, labels = coef_replace1, name = NULL, guide = "none") +
      facet_grid(trait ~ parameter, scales = "free_y", switch = "y",
                 labeller = labeller(trait = str_add_space, parameter = coef_replace1, .multi_line = FALSE)) +
      # facet_wrap(~ trait + parameter, scales = "free_y", nrow = length(traits), 
      #            labeller = labeller(trait = function(x) abbreviate(x, 2), parameter = coef_replace1, .multi_line = FALSE)) +
      theme_presentation(base_size = 10) +
      theme(strip.placement = "outside", legend.position = "bottom")
    
    return(g)
    
  })


# Prepare the plot list for combining
stability_mae_plotlist1 <- stability_mae_plotlist %>%
  map(~. + theme(axis.title = element_blank())) %>%
  map_at(-1:-3, ~ . + theme(strip.text.x = element_blank())) %>%
  map_at(-seq(1, length(traits) * 3, by = 3), ~. + theme(strip.text.y = element_blank())) %>%
  map_at(-seq(length(stability_mae_plotlist)-2, length(stability_mae_plotlist)), ~. + theme(axis.text.x = element_blank()))

# Merge, add axis grobs and legend
y_axis_grob1 <- grid::textGrob('Mean absolute deviation', gp = grid::gpar(fontsize = 10), rot = 90)
x_axis_grob <- grid::textGrob("Proportion of total environments", gp = grid::gpar(fontsize = 10))


main_plots <- wrap_plots(stability_mae_plotlist1, ncol = 3)
main_plots1 <- wrap_plots(main_plots, x_axis_grob, ncol = 1, heights = c(1, 0.02))
main_plots2 <- wrap_plots(y_axis_grob1, main_plots1, nrow = 1, widths = c(0.02, 1))


## Plot marker discovery rate
univariate_gwas_sigmar_df <- univariate_gwas_sigmar_grange %>%
  as.data.frame() %>%
  select(trait, parameter, marker)


# Factorize parameters
univariate_gwas_sample_sigmar_discovery_rate1 <- univariate_gwas_sigmar_df %>%
  left_join(., distinct(univariate_gwas_sample_discovery_rate, trait, pEnv)) %>%
  left_join(., univariate_gwas_sample_discovery_rate) %>%
  mutate(orig_signif_marker = TRUE,
         discovery_rate = ifelse(is.na(discovery_rate), 0, discovery_rate),
         parameter = factor(parameter, levels = names(coef_replace)))


# Create a plot of the discovery rate for QTL
# Only plot for originally significant markers
g_marker_discovery_rate <- univariate_gwas_sample_sigmar_discovery_rate1 %>%
  subset(parameter != "sd_d") %>%
  ggplot(aes(x = pEnv, y = discovery_rate, color = parameter)) +
  geom_line(aes(group = marker), linewidth = 0.5, alpha = 0.5) +
  geom_line(data = aggregate(discovery_rate ~ trait + parameter + pEnv, univariate_gwas_sample_sigmar_discovery_rate1, mean, subset = parameter != "sd_d"),
            linewidth = 1) +
  scale_color_manual(values = parameter_colors_use, guide = "none") +
  scale_y_continuous(name = "Marker-trait association discovery rate", breaks = pretty) +
  scale_x_continuous(name = "Proportion of total environments", limits = c(0.1, 0.9), breaks = seq(0.2, 0.8, by = 0.2)) +
  facet_grid(trait ~ parameter, labeller = labeller(trait = str_add_space, parameter = coef_replace1), switch = "y") +
  theme_presentation(base_size = 10) +
  theme(axis.title.x = element_blank())


# Combine plots
g_sampling_stability_mae_combined <- plot_grid(main_plots, g_marker_discovery_rate, 
                                               nrow = 1, labels = c("A", "B"), align = "v", axis = "tb", rel_widths = c(1, 0.6))
g_sampling_stability_mae_combined1 <- plot_grid(g_sampling_stability_mae_combined, x_axis_grob, ncol = 1, rel_heights = c(1, 0.05))

# Save
ggsave(filename = "figure3_sampling_results.jpg", plot = g_sampling_stability_mae_combined1, 
       path = fig_dir, width = 8, height = 6, dpi = 500)



# Example markers
univariate_gwas_sample_sigmar_discovery_rate1 %>%
  filter(pEnv < 0.5, discovery_rate > 0.5)

# Example markers
univariate_gwas_sample_sigmar_discovery_rate1 %>%
  filter(discovery_rate > 0.5, parameter == "b")


# Figure 4: genomic prediction accuracy -----------------------------------

# Part A: k-fold prediction accuracy
g_kfold_acc <- genotype_stability_kfold_acc %>%
  mutate(parameter = factor(parameter, levels = names(coef_replace))) %>%
  ggplot(aes(x = trait, y = accuracy, fill = parameter)) +
  geom_boxplot(alpha = 0.8, linewidth = 0.2, key_glyph = "rect") +
  scale_x_discrete(labels = as_mapper(~str_replace(str_add_space(.), " ", "\n")), name = NULL) +
  scale_y_continuous(name = "Cross-validation prediction accuracy", breaks = pretty) +
  scale_fill_manual(values = parameter_colors_use, labels = coef_replace1, name = NULL) +
  theme_presentation(base_size = 10) +
  theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        legend.position = "top")


# Part B - predicted vs observed values in POV
genotype_stability_pov1 <- genotype_stability_pov_results %>%
  unnest(results) %>%
  unnest(results)

# Calculate accuracy for annotations
genotype_stability_pov_acc_ann <- genotype_stability_pov_acc %>%
  mutate(annotation = format_numbers(accuracy, 2),
         annotation = paste0("italic(r)==", annotation),
         parameter = factor(parameter, levels = names(coef_replace)))

# List of plots
genotype_stability_pov_plotlist <- genotype_stability_pov1 %>%
  mutate(parameter = factor(parameter, levels = names(coef_replace))) %>%
  group_by(trait, parameter) %>%
  do(plot = {
    df <- .
    ggplot(df, aes(x = predicted.value, y = value, color = parameter)) +
      geom_point() +
      facet_grid(trait ~ ., labeller = labeller(trait = str_add_space), switch = "y") +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE) +
      geom_text(data = left_join(distinct(df, trait, parameter), genotype_stability_pov_acc_ann, by = c("trait", "parameter")),
                aes(x = Inf, y = -Inf, label = annotation), parse = TRUE, hjust = 1.2, vjust = -0.5, size = 2, inherit.aes = FALSE) +
      scale_color_manual(values = parameter_colors_use, guide = "none") +
      scale_x_continuous(name = "Predicted value", breaks = pretty, guide = guide_axis(check.overlap = TRUE)) +
      scale_y_continuous(name = "Observed value", breaks = pretty) +
      theme_presentation(base_size = 10) +
      theme(strip.placement = "outside")
  }) %>% ungroup()
  
# Modify the plotlist
genotype_stability_pov_plotlist1 <- genotype_stability_pov_plotlist %>%
  mutate(plot1 = map(.x = plot, ~ . + theme(axis.title = element_blank())),
         plot1 = map_if(.x = plot1, .p = parameter != "g", ~ . + theme(strip.text = element_blank())))

# Combine plots
g_genotype_stability_pov <- plot_grid(plotlist = genotype_stability_pov_plotlist1$plot1, nrow = length(traits), 
                                      rel_widths = c(1, 0.9, 0.9))

# Add y and x axis grobs
y_axis_grob1 <- grid::textGrob("Observed value", gp = grid::gpar(fontsize = 10), rot = 90)
x_axis_grob1 <- grid::textGrob("Predicted value", gp = grid::gpar(fontsize = 10), rot = 0)

g_genotype_stability_pov1 <- plot_grid(y_axis_grob1, plot_grid(g_genotype_stability_pov, x_axis_grob1, ncol = 1, rel_heights = c(1, 0.04)),
                                       nrow = 1, rel_widths = c(0.04, 1))





## Part C - prediction accuracy with environmental samples
# Plot
a <- 0.10

g_stability_samples_pov_acc <- genotype_stability_samples_pov_acc %>%
  group_by(trait, parameter, pEnv, target) %>%
  summarize_at(vars(accuracy), list(mean = mean, upper = ~quantile(., 1 - a/2, na.rm = T), lower = ~quantile(., a/2, na.rm = T)), 
               na.rm = TRUE) %>%
  ungroup() %>%
  filter(target == "predictions_test") %>%
  mutate(parameter = factor(parameter, levels = names(coef_replace))) %>%
  ggplot(aes(x = pEnv, y = mean, color = parameter)) +
  geom_hline(data = genotype_stability_pov_acc_ann, aes(yintercept = accuracy)) +
  geom_label(data = subset(genotype_stability_pov_acc_ann, trait == "GrainProtein" & parameter == "g"),
            aes(label = "Accuracy using all\nenvironments", x = 0.1, y = 0), inherit.aes = FALSE, hjust = 0, size = 2, label.size = 0) +
  geom_segment(data = subset(genotype_stability_pov_acc_ann, trait == "GrainProtein" & parameter == "g"),
               aes(x = 0.1, y = 0.1, xend = 0.1, yend = 0.5), inherit.aes = FALSE) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = parameter), alpha = 0.2, lwd = 0.5) +
  geom_line(lwd = 0.8) +
  # geom_point() +
  scale_x_continuous(name = "Proportion of total environments", breaks = pretty) +
  scale_y_continuous(name = "Between-generation prediction accuracy", breaks = pretty) +
  scale_fill_manual(values = parameter_colors_use, guide = "none") +
  scale_color_manual(values = parameter_colors_use, guide = "none") +
  facet_grid(trait ~ parameter, switch = "y", labeller = labeller(trait = str_add_space)) +
  theme_presentation(base_size = 10) +
  theme(strip.placement = "outside", strip.text.x = element_blank(), legend.position = "bottom")


# Combine the plots
left <- plot_grid(g_kfold_acc, g_stability_samples_pov_acc, ncol = 1, rel_heights = c(0.7, 1), labels = c("A", "C"))
right <- g_genotype_stability_pov1

g_combined <- plot_grid(left, right, nrow = 1, rel_widths = c(0.65, 1), labels = c("", "B"))


# Save
ggsave(filename = "figure4_genomewide_prediction_results.jpg", plot = g_combined, 
       path = fig_dir, width = 10, height = 7, dpi = 500)




# Find the proportion of environments in which prediction accuracy is 80% or more of
# the accuracy when using all data
genotype_stability_samples_pov_acc %>%
  group_by(trait, parameter, pEnv, target) %>%
  summarize(accuracy = mean(accuracy, na.rm = T), .groups = "drop") %>%
  left_join(., select(genotype_stability_pov_acc, trait, parameter, full_acc = accuracy)) %>%
  mutate(p80orMore = (accuracy / full_acc) >= 0.80) %>%
  filter(target == "predictions_test", p80orMore) %>%
  arrange(trait, parameter, pEnv) %>%
  group_by(trait, parameter) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  arrange(parameter, pEnv)







# Table 1: Stability summary ----------------------------------------------


# Create a table of summary statistics for stability, including heritability and
# correlations

# Heritability table
genotype_stability_herit_toprint <- genotype_stability_herit %>%
  mutate_at(vars(h2, std_error), ~ format_numbers(., signif.digits = 2)) %>%
  mutate(parameter = factor(parameter, levels = names(coef_replace)),
         annot = paste0(h2, " (", std_error, ")")) %>%
  select(trait, parameter, annot) %>%
  spread(parameter, annot)

# Genetic correlation tables
genotype_stability_cor_toprint <- genotype_stability_cor %>%
  select(trait, var_pair, correlation) %>% 
  unnest(correlation) %>%
  filter(parameter == "corG", is.finite(estimate)) %>%
  mutate_at(vars(estimate, std_error), ~ format_numbers(., signif.digits = 2)) %>%
  mutate(annot = paste0(estimate, " (", std_error, ")")) %>%
  select(trait, var_pair, annot) %>%
  spread(var_pair, annot)

# Combine
genotype_stability_herit_cor_toprint <- left_join(genotype_stability_herit_toprint, genotype_stability_cor_toprint) %>%
  mutate(trait = str_add_space(trait))
  
# Save
write_csv(x = genotype_stability_herit_cor_toprint, file = file.path(fig_dir, "table1_stability_heritability_correlation.csv"))
  



# Table 2: Summary of marker-trait associations ---------------------------

univariate_gwas_sigmar_df <- as.data.frame(univariate_gwas_sigmar_grange) %>%
  left_join(univariate_gwas_farm_bestModel_qq, .)

# Summarize the number of PCs and associations per trait
univariate_gwas_sigmar_summ1 <- univariate_gwas_sigmar_df %>%
  group_by(trait, parameter) %>%
  summarize(nPC = unique(bestModel_qq_nPC), nMTA = sum(!is.na(marker)), .groups = "drop")
  

## Calculate the min/max/mean distance of each MTA with T3 qtl

# First update t3_granges to look at relevant traits
t3_qtl_granges1 <- t3_qtl_granges %>%
  as.data.frame() %>%
  subset(t3_trait %in% c("grain protein", "grain yield", "heading date", "heading date (Julian)", "plant height", "test weight")) %>%
  mutate(t3_trait1 = case_when(
    t3_trait %in% c("diastatic power", "grain protein", "alpha amylase") ~ "GrainProtein",
    t3_trait %in% c("grain yield") ~ "GrainYield",
    t3_trait %in% c("heading date", "heading date (Julian)") ~ "HeadingDate",
    t3_trait %in% c("plant height") ~ "PlantHeight",
    t3_trait %in% c("plump grain", "breeders plump grain") ~ "PlumpGrain",
    t3_trait %in% c("test weight") ~ "TestWeight"
  )) %>%
  GRanges()

# Now find the distance from each MTA to the nearest QTL
univariate_gwas_sigmar_nearestQTL <- list()

for (tr in traits) {
  x <- subset(univariate_gwas_sigmar_grange, trait == tr)
  subject <- subset(t3_qtl_granges1, t3_trait1 == tr)
  
  distNearest <- distanceToNearest(x = x, subject = subject) %>%
    subset(distance <= sig_snp_window)
    # Pull out the nearest QTL
    
  x1 <- cbind(as.data.frame(x[from(distNearest)]), distance_to_nearest_qtl = as.data.frame(distNearest)$distance)
  univariate_gwas_sigmar_nearestQTL[[tr]] <- x1
  
}

# Calculate the mean, min, max distance of each significant SNP to the nearest T3 QTL
univariate_gwas_overlap_qtl_distance <- bind_rows(univariate_gwas_sigmar_nearestQTL) %>%
  group_by(trait, parameter) %>%
  summarize_at(vars(distance_to_nearest_qtl), list(median_dist = median, min_dist = min, max_dist = max)) %>%
  ungroup()

univariate_gwas_T3Overlap_count <- mergeByOverlaps(query = univariate_gwas_sigmar_grange, subject = t3_qtl_granges1, maxgap = sig_snp_window) %>% 
  subset(trait == t3_trait1)
univariate_gwas_T3Overlap_count <- as.data.frame(unique(univariate_gwas_T3Overlap_count[,c("marker", "trait", "parameter")])) %>%
  mutate(MTA_T3_overlap = 1)


# Count the number of MTA with overlap with nearby T3 QTL
univariate_gwas_T3overlap_summ2 <- as.data.frame(univariate_gwas_sigmar_grange) %>%
  left_join(., univariate_gwas_T3Overlap_count) %>%
  mutate(MTA_T3_overlap = ifelse(is.na(MTA_T3_overlap), 0, MTA_T3_overlap)) %>%
  group_by(trait, parameter) %>%
  summarize(nMTA_T3_overlap = sum(MTA_T3_overlap), .groups = "drop") %>%
  left_join(univariate_gwas_sigmar_summ1, .) %>%
  mutate(nMTA_T3_overlap = ifelse(is.na(nMTA_T3_overlap), 0, nMTA_T3_overlap))


# Now find the distance from each MTA to the nearest QTL
distNearestGene <- distanceToNearest(x = univariate_gwas_sigmar_grange, subject = unique(univariate_gwas_sigmar_grange_geneOverlap$barley_txdb2_genes_grange)) %>%
  subset(distance <= sig_snp_window)

univariate_gwas_sigmar_nearestGene <- univariate_gwas_sigmar_grange[from(distNearestGene)] %>%
  as.data.frame() %>%
  mutate(distance_to_nearest_gene =  as.data.frame(distNearestGene)$distance,
         nearest_gene_id = unique(univariate_gwas_sigmar_grange_geneOverlap$barley_txdb2_genes_grange)[to(distNearestGene)]$gene_id,
         nearest_gene_name = unique(univariate_gwas_sigmar_grange_geneOverlap$barley_txdb2_genes_grange)[to(distNearestGene)]$Name)
  

# Calculate the mean, min, max distance of each significant SNP to the nearest gene
univariate_gwas_overlap_gene_distance <- univariate_gwas_sigmar_nearestGene %>%
  group_by(trait, parameter) %>%
  summarize_at(vars(distance_to_nearest_gene), list(median_dist = median, min_dist = min, max_dist = max)) %>%
  ungroup()


# Same thing for genes
univariate_gwas_geneOverlap_summ1 <- univariate_gwas_sigmar_grange_geneOverlap %>%
  as.data.frame() %>%
  group_by(trait, parameter) %>%
  summarize(nMTA_geneOverlap = n_distinct(marker), .groups = "drop") %>%
  left_join(univariate_gwas_sigmar_summ1, .) %>%
  mutate(nMTA_geneOverlap = ifelse(is.na(nMTA_geneOverlap), 0, nMTA_geneOverlap))


# Table of variance explained
univariate_gwas_varexp <- univariate_gwas_farm_varexp %>%
  filter(kinship == "farmcpu") %>%
  left_join(select(univariate_gwas_farm_bestModel_qq, trait, parameter, nPC = bestModel_qq_nPC), .) %>%
  unnest(r2_sigmar_df) %>%
  group_by(trait, parameter, nPC) %>%
  mutate_at(vars(r2), list(r2_mean = mean, r2_min = min, r2_max = max)) %>%
  ungroup() %>%
  mutate_at(vars(contains("r2")), format_numbers, 2) %>%
  select(-marker, -r2) %>%
  distinct() %>%
  mutate(r2_sigmar_range = paste0(r2_min, "-", r2_max))

# Table of variance explained for each SNP
univariate_gwas_varexp_persnp <- univariate_gwas_farm_varexp %>%
  filter(kinship == "farmcpu") %>%
  left_join(select(univariate_gwas_farm_bestModel_qq, trait, parameter, nPC = bestModel_qq_nPC), .) %>%
  unnest(r2_sigmar_df) %>%
  select(trait, parameter, marker, r2)



# Combine
table2_toprint <- univariate_gwas_sigmar_summ1 %>%
  left_join(., select(univariate_gwas_varexp, trait, parameter, nPC, r2_kinship, r2_sigmar, r2_sigmar_range)) %>%
  # Add overlap with T3 QTL
  left_join(., univariate_gwas_T3overlap_summ2) %>%
  # # Add distance to QTL
  # left_join(., univariate_gwas_overlap_qtl_distance) %>%
  # mutate_at(vars(contains("dist")), ~format_numbers(. / 1000, 2) %>% str_remove(string = ., pattern = "\\.$")) %>%
  # # Clean up
  # mutate(nearest_qtl_kb = paste0(median_dist, " (", min_dist, ", ", max_dist, ")"),
  #        nearest_qtl_kb = ifelse(grepl(pattern = "NA", x = nearest_qtl_kb), as.character(NA), nearest_qtl_kb)) %>%
  # select(-median_dist, -min_dist, -max_dist) %>%
  # Add gene overlap
  left_join(., univariate_gwas_geneOverlap_summ1) %>%
  # # Add gene distance
  # left_join(., univariate_gwas_overlap_gene_distance) %>%
  # mutate_at(vars(contains("dist")), ~format_numbers(. / 1000, 2) %>% str_remove(string = ., pattern = "\\.$")) %>%
  # # Clean up
  # mutate(nearest_gene_kb = paste0(median_dist, " (", min_dist, ", ", max_dist, ")"),
  #        nearest_gene_kb = ifelse(grepl(pattern = "NA", x = nearest_gene_kb), as.character(NA), nearest_gene_kb)) %>%
  # select(-median_dist, -min_dist, -max_dist) %>%
  select(names(.))



# Save table 2
table2_toprint %>%
  mutate(parameter = factor(parameter, levels = names(coef_replace))) %>%
  arrange(trait, parameter) %>%
  write_csv(x = ., file = file.path(fig_dir, "table2_gwas_summary.csv"))






# Supplemental Table S1: All significant marker-trait associations --------

# Also display distance to nearest gene/QTL and the information for that gene/QTL
all_gwas_sigmar <- univariate_gwas_sigmar_grange %>%
  as.data.frame() %>%
  select(trait, parameter, marker, chrom = seqnames, position = start, p_value, maf, effect)

t3_qtl_granges2 <- subset(t3_qtl_granges1, phenotyping_source != "experiment_set")
sigmar_qtl_nearest_distance <- distanceToNearest(univariate_gwas_sigmar_grange, t3_qtl_granges2)

all_gwas_sigmar_nearest_qtl <- t3_qtl_granges2[to(sigmar_qtl_nearest_distance)]
all_gwas_sigmar_nearest_qtl$marker <- univariate_gwas_sigmar_grange[from(sigmar_qtl_nearest_distance)]$marker
all_gwas_sigmar_nearest_qtl$trait <- univariate_gwas_sigmar_grange[from(sigmar_qtl_nearest_distance)]$trait
all_gwas_sigmar_nearest_qtl$parameter <- univariate_gwas_sigmar_grange[from(sigmar_qtl_nearest_distance)]$parameter
all_gwas_sigmar_nearest_qtl$distance <- sigmar_qtl_nearest_distance@elementMetadata$distance

nOverlapsPerMarker <- findOverlaps(query = univariate_gwas_sigmar_grange, subject = t3_qtl_granges2, maxgap = sig_snp_window) %>% 
  as.data.frame() %>% 
  split(.$queryHits) %>% 
  map_dbl(nrow)


# Data frame of number of nearby QTL and nearest QTL
all_gwas_sigmar1 <- all_gwas_sigmar_nearest_qtl %>%
  as.data.frame(row.names = NULL) %>%
  select(trait, parameter, marker, nearest_qtl_pos = start, nearest_qtl_distance = distance) %>%
  filter(nearest_qtl_distance <= sig_snp_window) %>%
  left_join(all_gwas_sigmar, .) %>%
  mutate(nOverlappingQTL = nOverlapsPerMarker)


# Do the same thing for nearby genes
barley_txdb2_genes_grange <- univariate_gwas_sigmar_grange_geneOverlap$barley_txdb2_genes_grange
sigmar_gene_nearest_distance <- distanceToNearest(univariate_gwas_sigmar_grange, barley_txdb2_genes_grange)

all_gwas_sigmar_nearest_gene <- barley_txdb2_genes_grange[to(sigmar_gene_nearest_distance)]
all_gwas_sigmar_nearest_gene$marker <- univariate_gwas_sigmar_grange[from(sigmar_gene_nearest_distance)]$marker
all_gwas_sigmar_nearest_gene$trait <- univariate_gwas_sigmar_grange[from(sigmar_gene_nearest_distance)]$trait
all_gwas_sigmar_nearest_gene$parameter <- univariate_gwas_sigmar_grange[from(sigmar_gene_nearest_distance)]$parameter
all_gwas_sigmar_nearest_gene$distance <- sigmar_gene_nearest_distance@elementMetadata$distance

nGeneOverlapsPerMarker <- findOverlaps(query = univariate_gwas_sigmar_grange, subject = barley_txdb2_genes_grange, maxgap = sig_snp_window) %>% 
  as.data.frame() %>% 
  split(.$queryHits) %>% 
  map_dbl(nrow)

# Data frame of number of nearby QTL and nearest QTL
all_gwas_sigmar2 <- all_gwas_sigmar_nearest_gene %>%
  as.data.frame(row.names = NULL) %>%
  select(trait, parameter, marker, nearest_gene_pos_start = start, nearest_gene_pos_end = end, gene_id, gene_name = Name, 
         nearest_gene_distance = distance) %>%
  filter(nearest_gene_distance <= sig_snp_window) %>%
  left_join(all_gwas_sigmar1, .) %>%
  mutate(nOverlappingGenes = nGeneOverlapsPerMarker) %>%
  left_join(., univariate_gwas_varexp_persnp) %>%
  select(trait, parameter, marker, chrom, position, p_value, maf, effect, r2, names(.)) %>%
  # Add named genes
  mutate(.id = seq_len(nrow(.)),
         chrom1 = sub("H", "", chrom)) %>%
  split(.$.id) %>%
  map_df(~{
    genes1 <- subset(barley_known_genes, chrom1 == .x$chrom1 & start >= .x$position - sig_snp_window & start <= .x$position + sig_snp_window)
    genes2 <- subset(gene_overlap_positions, chrom1 == .x$chrom1 & pos >= .x$position - sig_snp_window & pos <= .x$position + sig_snp_window)
    
    # Add these genes and their distances
    genes1 <- genes1 %>%
      select(chrom1, nearest_named_gene_pos_start = start, nearest_named_gene_pos_end = end, nearest_named_gene_name = gene) %>%
      mutate(nearest_named_gene_distance = pmin(abs(nearest_named_gene_pos_start - .x$position), abs(nearest_named_gene_pos_end - .x$position))) %>%
      distinct() %>%
      top_n(n = 1, wt = -nearest_named_gene_distance)
    genes2 <- genes2 %>%
      distinct() %>%
      select(chrom1, nearest_named_gene_pos_start = start, nearest_named_gene_pos_end = end, nearest_named_gene_name = Name) %>%
      mutate(nearest_named_gene_distance = pmin(abs(nearest_named_gene_pos_start - .x$position), abs(nearest_named_gene_pos_end - .x$position))) %>%
      top_n(n = 1, wt = -nearest_named_gene_distance)
    
    if (nrow(genes1) > 0 & nrow(genes2) > 0) {
      if (genes1$nearest_named_gene_distance < genes2$nearest_named_gene_distance) {
        merge(.x, genes1)
      } else {
        merge(.x, genes2)
      }
    } else if (nrow(genes1) > 0) {
      merge(.x, genes1)
    } else if (nrow(genes2) > 0) {
      merge(.x, genes2)
    } else {
      .x
    }
  })



# Look up GO terms for nearby genes ---------------------------------------

library(GO.db)
library(clusterProfiler)

# Load gene annotations
gene_annotations <- read_tsv(list.files(path = data_dir, pattern = "annotation", full.names = TRUE)) %>%
  rename_all(~tolower(make.names(.))) %>%
  select(-...8)

## Get the names and GO terms for all genes within the set range of each marker
gene_ids_annotations <- univariate_gwas_sigmar_grange_geneOverlap %>%
  as.data.frame() %>%
  distinct(gene_id) %>% 
  left_join(., gene_annotations)

# Look up the annotations for the closest genes
all_gwas_sigmar3 <- all_gwas_sigmar2 %>%
  left_join(., gene_ids_annotations)


# Get unique go terms
go_terms <- unique(unlist(strsplit(all_gwas_sigmar3$go.ids, ", ")))
go_terms <- subset(go_terms, go_terms != "none")

# Get the descriptions for those GO terms
all_go <- as.list(GOTERM)

# Build a DF
go_term_descr <- NULL
for (go in go_terms) {
  if (go %in% names(all_go)) {
    term <- Term(all_go[[go]])
    def <- Definition(all_go[[go]])
    go_term_descr <- rbind(go_term_descr, data.frame(go_id = go, go_term = term, go_definition = def))
  }
}


# Save this table
tableS1 <- all_gwas_sigmar3 %>%
  mutate(parameter = coef_replace1[parameter]) %>%
  dplyr::select(marker, chrom, position, trait, parameter, p_value, maf, alt_allele_effect = effect, r2, n_overlapping_qtl = nOverlappingQTL,
         nearest_qtl_distance, n_overlapping_genes = nOverlappingGenes, nearest_named_gene_distance, nearest_named_gene_name,  
         nearest_gene_distance, nearest_gene_id = gene_id, nearest_gene_description = description, nearest_gene_go_id = go.ids)

tableS1a <- tableS1 %>%
  # Split open go ids
  as_tibble() %>%
  mutate(nearest_gene_go_id = strsplit(nearest_gene_go_id, ", ")) %>%
  unnest(nearest_gene_go_id) %>%
  left_join(., go_term_descr, by = c("nearest_gene_go_id" = "go_id"))



# Save
write_csv(x = tableS1a, file = file.path(fig_dir, "tableS1_marker_trait_associations.csv"), na = "")


# Cluster go terms
go_data <- godata(ont = "BP")
go_terms <- unique(tableS1a$nearest_gene_go_id)
go_terms <- subset(go_terms, go_terms != "none")
sim_matrix <- mgoSim(go_terms, go_terms, semData = go_data, measure = "Wang", combine = NULL)
hc <- hclust(as.dist(1 - sim_matrix))  # Convert similarity to distance
plot(hc)
groups <- cutree(tree = hc, h = 0.5)

# Run a quick comparison of go terms for the three different parameters
tableS1a %>%
  group_by(parameter, nearest_gene_go_id, go_term) %>%
  count() %>%
  ungroup() %>%
  arrange(go_term, parameter) %>%
  view()

# There does not seem to be any significant enrichment for certain go terms between the 
# mean and plasticity results



# Other results referenced in the paper -----------------------------------

intersectionPoint <- function(l1, l2){
  x <- (l2[1] - l1[1]) / (l1[2] - l2[2])
  y <- l1[1] + l1[2] * x
  return(xy=c(x, y))
}

intersectInRange <- function(l1, l2, xlim) {
  xy <- intersectionPoint(l1 = l1, l2 = l2)
  xlim[1] <= xy[1] & xy[1] <= xlim[2]
}


# Calculate the number of pairwise crossovers

genotype_stability_crossovers <- genotype_stability %>%
  select(-ranova) %>%
  left_join(., nest(group_by(genotype_environment_mean, trait))) %>%
  group_by(trait) %>%
  do({
    row <- .
    df <- row$regression[[1]] %>%
      filter(line_name %in% tp)
    
    # Pairwise combination of lines 
    crossover <- combn(x = as.character(df$line_name), m = 2, FUN = function(genos) {
      df1 <- subset(df, line_name %in% genos)
      intersectInRange(l1 = c(df1$g[1], df1$b[1]), l2 = c(df1$g[2], df$b[2]), xlim = range(row$data[[1]]$h))
    })
    
    tibble(pCrossovers = mean(crossover))

  }) %>% ungroup()


genotype_stability_crossovers %>%
  filter(trait %in% traits)

# # A tibble: 5 x 2
# trait        pCrossovers
# 1 GrainProtein       0.678
# 2 GrainYield         0.444
# 3 HeadingDate        0.201
# 4 PlantHeight        0.340
# 5 TestWeight         0.375




# Calculate the percentage of MTA for genotype mean / stability that overlap

univariate_gwas_sigmar_grange_df <- univariate_gwas_sigmar_grange %>%
  as.data.frame() %>%
  group_by(trait, parameter) %>%
  nest() %>%
  ungroup()

univariate_gwas_sigmar_grange_df %>%
  crossing(., ., .name_repair = tidyr_legacy) %>%
  filter(trait == trait1, parameter != parameter1) %>%
  select(-trait1) %>%
  group_by(trait, parameter, parameter1) %>%
  do({
    row <- .
    gr1 <- GRanges(row$data[[1]])
    gr2 <- GRanges(row$data1[[1]])
    
    # Find the overlaps
    overlaps <- findOverlaps(query = gr1, subject = gr2, maxgap = sig_snp_window)
    # Find exact overlaps
    exactOverlaps <- findOverlaps(query = gr1, subject = gr2, maxgap = 0)
    
    # Return
    tibble(nMar1 = length(gr1), nMar2 = length(gr2), nOverlaps = length(overlaps), nExactOverlaps = length(exactOverlaps))
    
  })


## Plot LD
map <- snp_info_all %>% 
  subset(marker %in% colnames(geno_mat_LD), c("marker", "chrom", "pos"))

ld_plot(geno = geno_mat_LD, map = map)




