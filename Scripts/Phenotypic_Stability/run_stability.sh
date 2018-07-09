#!/bin/bash

#PBS -l walltime=12:00:00,mem=24gb,nodes=1:ppn=16
#PBS -N stability_genetic_correlation_permutation
#PBS -M neyha001@umn.edu
#PBS -m abe
#PBS -r n

# Change the working directory
cd /panfs/roc/groups/6/smithkp/neyha001/QTLMapping/S2MET_Mapping/Scripts/Phenotypic_Stability

module load R/3.2.0_intel_mkl

# Marker effect stability variance component test
# Rscript S2MET_marker_stability_varcomp.R

# Marker effect by environment computation
# Rscript S2MET_marker_stability_computation.R

# # Stability prediction
# Rscript stability_prediction.R

# Stability genetic correlation
Rscript stability_analysis_permutation.R
