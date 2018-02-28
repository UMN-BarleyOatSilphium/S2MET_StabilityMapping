## Download and compile GWAS data from the Triticeae Toolbox
## 
## This script will outline the steps of downloading GWAS results from the Triticeae
## Toolbox database (T3), compiling the marker data, and comparing the results
## of our analysis to these previous GWAS results
## 

# Load some packages
library(tidyverse)
library(readxl)
library(GenomicRanges)

# Directory containing the project repository
repo_dir <- getwd()

# Project and other directories
source(file.path(repo_dir, "source.R"))

# Read in the QTL metadata files that were downloaded from T3
qtl_meta <- map_df(list.files(data_dir, pattern = "meta", full.names = TRUE), read_csv)

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

# Read in the significant markers identified in our analysis
load(file.path(result_dir, "gwas_mlmm_marker_associations.RData"))

## Create GRanges for each marker data.frame
qtl_meta_grange <- qtl_meta_use %>%
  split(.$trait) %>%
  map(~makeGRangesFromDataFrame(., keep.extra.columns = TRUE, start.field = "position", 
                                end.field = "position"))





# Add a window of 5 Mbp to either side of the significant markers
window <- 5000000

gwas_mlmm_grange <- gwas_mlmm_marker %>%
  mutate(start = pos - window, end = pos + window) %>%
  split(.$trait) %>%
  map(~makeGRangesFromDataFrame(., keep.extra.columns = TRUE))

## For each trait, see how many of the markers detected in our analysis overlap
## with those on T3
gwas_markers_overlap <- list(gwas_mlmm_grange, qtl_meta_grange) %>%
  pmap_df(function(gwas, t3) {
    overlaps <- findOverlaps(query = gwas, subject = t3)
    
    # Add the overlapping markers to the gwas marker grange
    overlap_data <- bind_cols(
      as.data.frame(gwas[queryHits(overlaps)]),
      as.data.frame(t3[subjectHits(overlaps)])
    )
    
    # Bind with the markers that did not overlaps
    non_overlap_data <- as.data.frame(gwas[setdiff(seq_along(gwas), queryHits(overlaps))])
    
    # Sort and remove superfluous columns
    bind_rows(overlap_data, non_overlap_data) %>%
      select(trait, coef, marker, chrom = seqnames, pos, start, end, MAF:WA, 
             qtl_marker = marker1, qtl_pos = start1, gene, feature) %>% 
      arrange(trait, coef, chrom, pos) %>%

    
  })

## Write to a table
save_file <- file.path(fig_dir, "gwas_significant_associations.csv")
write_csv(x = gwas_markers_overlap, path = save_file)


# Read in the annotated excel file
gwas_markers_overlap <- read_excel(path = file.path(fig_dir, "gwas_significant_associations.xlsx"),
                                   na = c("NA")) %>%
  mutate(qtl_pos = parse_number(qtl_pos))

# Number of associations per trait
gwas_markers_overlap


## Count the number of overlaps per significant GWAS marker
gwas_markers_overlap_count <- gwas_markers_overlap %>% 
  group_by(trait, coef, marker) %>% 
  summarize(n_overlap = sum(!is.na(source)))

## Summarize overlaps per trait-coef combination
gwas_markers_overlap_count %>% 
  summarize(count_overlap = sum(n_overlap > 0),
            prop_overlap = mean(n_overlap > 0))

# Average distance of overlapping markers with the closest overlapping QTL
gwas_markers_overlap %>% 
  # Calculate the distance from the GWAS marker to the QTL marker
  mutate(marker_qtl_distance = abs(pos - qtl_pos)) %>%
  group_by(trait, coef, marker) %>% 
  summarize(marker_qtl_distance = min(marker_qtl_distance)) %>%
  summarize(mean_distance = mean(marker_qtl_distance, na.rm = TRUE))












