<<<<<<< HEAD
# Scripts
 
The order and description of the scripts is as follows:

`01_calculate_mean_and_stability.R` - calculates stability using all data or subsets of data.

`02_stability_mapping.R` - run genomewide association analysis for the genotype mean and stability estimates.

`03_stability_mapping_resampling.R` - run genomewide association analysis using the subsets of environments.

`04_stability_mapping_analysis.R` - analysis of mapping results, including identification of overlapping QTL and genes.

`05_stability_prediction.R` - run genomewide prediction for the genotype mean and stability.

`06_figures.R` - create figures.
=======
# Scripts
 
This script includes information on the order in which scripts in this project should be executed.
It also will check the datestamp of RData files to see check what analyses have been run recently

## Order of scripts
Below, R code is provided that will run scripts in their entirety. The order of the scripts is the
recommended order for re-running the analysis

## Part 1 - Phenotypic data analysis
source(file.path(script_dir, "phenotype_data_summary.R"))   # Summarizes phenotype data across environments

## Part 2 - Calculating genotype mean and stability
source(file.path(script_dir, "PhenotypicStability/trait_mean_stability_computation.R"))  # Calculates trait mean and stability estimates
source(file.path(script_dir, "PhenotypicStability/marker_effect_by_env_and_stability.R"))    # Calculates marker effect stability
source(file.path(script_dir, "PhenotypicStability/stability_marker_varcomp.R")) # Calculates the proportion of variance in mean and stability attributed to genomewide markers or subsets

## Part 3 - Association Mapping
source(file.path(script_dir, "AssociationMapping/stability_mapping.R"))   # Main script for association mapping and MLMM
source(file.path(script_dir, "AssociationMapping/stability_mapping_resampling.R"))    # Association mapping using environmental subsets (required HPC resources)
source(file.path(script_dir, "AssociationMapping/pleiotropy_mapping.R"))    # Pleiotropy association mapping

## Part 4 - Genomewide predictions
source(file.path(script_dir, "GenomewidePrediction/stability_prediction.R"))    # Runs cross-validation and parent-offspring validation




## Analysis
## These scripts generally correspond to the scripts above, but I recommend running them after all of the "experiment"
## scripts are completed. Many of these scripts also generate figures
# source(file.path(script_dir, "PhenotypicStability/stability_analysis.R"))   # Analyze the stability results
# source(file.path(script_dir, "PhenotypicStability/marker_stability_analysis.R"))    # Analyze the marker effect stability results

# source(file.path(script_dir, "AssociationMapping/stability_mapping_analysis.R"))    # Analyze all mapping results

# source(file.path(script_dir, "GenomewidePrediction/stability_prediction_analysis.R"))   # Analyze the prediction results




>>>>>>> c042eb08f49fab9587e31a12e4f60ec19baa31ee
