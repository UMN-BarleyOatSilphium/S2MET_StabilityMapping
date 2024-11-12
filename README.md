
<!-- README.md is generated from README.Rmd. Please edit that file -->

# S2MET_StabilityMapping

## Description

This repository contains information and code for replicating the
analyses described in the paper:

Article Title: *Genomewide association and prediction of phenotypic
stability in barley*  
Journal: *The Plant Genome* (under review)  
Authors: Jeffrey L. Neyhart, Lucia Gutierrez, and Kevin P. Smith  
Article doi: TBD

## Navigation

### Data

Data used in this study are available from the [Triticeae
Toolbox](https://triticeaetoolbox.org/barley). See [this
README](git@github.com:neyhartj/S2MET_StabilityMapping/tree/master/Data)
for instructions on accessing this data.

### Code

Scripts used to complete the analysis and generate figures outlined in
the article above are available in the “Scripts” subfolder. See [this
README](git@github.com:neyhartj/S2MET_StabilityMapping/tree/master/Scripts)
for information on the scripts and their intended execution order.

Three scripts in this directory are used by all other scripts:

1.  `startub.R` - loads packages, creates directory links, and loads
    data.
2.  `functions.R` - loads additional functions into the environment.

The `06_figures.R` script in the “Scripts” subfolder produces the
figures found in the paper.

## Software/package versions

*R* version 4.4.1 was used to complete analyses in this project.

The following packages (with version information) were used to complete
this project:

| package         | version  |
|:----------------|:---------|
| FarmCPUpp       | 1.2.0    |
| GO.db           | 3.19.1   |
| GWASTools       | 1.50.0   |
| GenomicRanges   | 1.56.2   |
| biganalytics    | 1.1.22   |
| bigmemory       | 4.6.4    |
| blme            | 1.0-6    |
| clusterProfiler | 4.12.6   |
| cowplot         | 1.1.3    |
| lme4            | 1.1-35.5 |
| lmerTest        | 3.1-3    |
| modelr          | 0.1.11   |
| neyhart         | 0.1.1    |
| patchwork       | 1.3.0    |
| readxl          | 1.4.3    |
| rrBLUP          | 4.6.3    |
| scam            | 1.2-17   |
| snpStats        | 1.54.0   |
| sommer          | 4.3.6    |
| tidyverse       | 2.0.0    |
