## Use Large-scale T3 data
## 
## 
## 

library(lme4)
library(modelr)
library(rrBLUP)

# Source data
repo_dir <- getwd()
source(file.path(repo_dir, "source.R"))

## Read in the T3 data
T3_dir <- file.path(data_dir, "BCAP_Data")

# Phenotypes
pheno <- (list.files(T3_dir, pattern = ".txt") %>% .[!str_detect(., "hmp")]) %>%
  map_df(~read_tsv(file = file.path(T3_dir, .)) %>% gather(trait, value, -line, -trial)) %>%
  mutate(line = str_remove_all(line, "'"),
         trait = str_remove_all(str_to_title(trait), " ")) %>%
  filter(!is.na(value))


## First remove trials with < 50 entries
min_entries <- 50

pheno1 <- pheno %>% 
  group_by(trait, trial) %>%
  mutate(nObs = n_distinct(line)) %>%
  filter(nObs >= min_entries)

# Second filter lines that were observed in less than 5 environments
min_env <- 5

pheno2 <- pheno1 %>%
  group_by(line) %>%
  mutate(nObs = n_distinct(trial)) %>%
  filter(nObs >= min_env) %>%
  ungroup()

## Tidy for modeling
pheno_tomodel <- pheno2 %>%
  select(-nObs) %>%
  mutate_at(vars(line, trial), as.factor)



## Genotypes
geno <- read_tsv(file = file.path(T3_dir, "genotype.hmp.txt")) %>%
  filter(!chrom %in% c("Un", "UNK")) %>%
  mutate(chrom = parse_number(chrom))

# Convert to marker matrix
geno_mat <- t(select(geno, -rs:-pos))
geno_mat <- geno_mat[row.names(geno_mat) %in% levels(pheno_tomodel$line),]
colnames(geno_mat) <- geno$rs

# Missing markers
marker_miss <- colMeans(is.na(geno_mat)); hist(marker_miss)
geno_mat1 <- geno_mat[,marker_miss <= 0.80]

# Marker MAF
marker_af <- (colMeans(geno_mat1 + 1, na.rm = TRUE) / 2)
marker_maf <- pmin(marker_af, 1 - marker_af); hist(marker_maf)
geno_mat2 <- geno_mat1[,marker_maf >= round(15 / nrow(geno_mat1), 3)]

# Line missing
line_miss <- rowMeans(is.na(geno_mat2)); hist(line_miss)
geno_mat3 <- geno_mat2[line_miss <= 0.80, ]


## Impute
impute_out <- A.mat(X = geno_mat3, min.MAF = 0, max.missing = 1, return.imputed = TRUE)
geno_mat_impute <- impute_out$imputed

## Convert back to hmp
geno1 <- cbind(subset(geno, rs %in% colnames(geno_mat_impute), c(rs, chrom, pos)), t(geno_mat_impute))


## PCA
geno_pca <- prcomp(geno_mat_impute)





## Mixed model to calculate genotype and environment effect
pheno_fit1 <- pheno_tomodel %>%
  group_by(trait) %>%
  do(fit1 = lmer(value ~ (1|line) + trial, data = ., contrasts = list(trial = "contr.sum"))) %>%
  ungroup()

# Get the environmental effects
coefs <- fixef(fit1)
env_eff <- coefs[2:n_distinct(pheno_tomodel$trial)]

env_eff_df <- pheno_fit1 %>%
  mutate(coef = map(fit1, fixef),
         env_eff = map(coef, ~c(., -sum(.))),
         env_eff_df = map(env_eff, )
  
  
  
  data_frame(trial = levels(pheno_tomodel$trial), h = c(env_eff, -sum(env_eff)))

## FW regression
pheno_tomodel1 <- left_join(pheno_tomodel, env_eff_df)

fw_reg <- pheno_tomodel1 %>%
  rename(line_name = line, environment = trial) %>%
  group_by(trait, line_name) %>%
  do(calc_stability(., remove.outliers = FALSE)) %>%
  ungroup() %>%
  filter(type == "outliers") %>%
  mutate(g = map_dbl(model, ~coef(.)[1]))

## Convert to phenotypes for GWAS
fw_tomap <- fw_reg %>% 
  select(trait, line_name, b, g) %>%
  split(.$trait) %>%
  map(~select(., -trait))

## Map
gwas_out <- fw_tomap %>%
  map(~{
    df <- .
    geno_use <- geno1[,c(1:5, which(colnames(geno1) %in% df$line_name))]
    df1 <- filter(df, line_name %in% colnames(geno_use)) %>%
      as.data.frame()
    GWAS(pheno = df1, geno = geno_use, n.PC = 1, min.MAF = 0, P3D = TRUE, plot = FALSE)
  })
  
## QQ plots
library(qqman)

qq(pvector = gwas_out$GrainYield$g)
qq(pvector = gwas_out$GrainYield$b)

manhattan(x = gwas_out$GrainYield, chr = "chrom", bp = "pos", p = "g", snp = "rs")
manhattan(x = gwas_out$GrainYield, chr = "chrom", bp = "pos", p = "b", snp = "rs")
























