#!/bin/bash

#PBS -l walltime=12:00:00,mem=24gb,nodes=1:ppn=16
#PBS -N stability_genetic_correlation_permutation
#PBS -M neyha001@umn.edu
#PBS -m abe
#PBS -r n

# Change the working directory
cd /panfs/roc/groups/6/smithkp/neyha001/QTLMapping/S2MET_Mapping/Scripts/GenomewidePrediction

module load R/3.2.0_intel_mkl

# # Stability prediction
# Rscript stability_prediction.R
