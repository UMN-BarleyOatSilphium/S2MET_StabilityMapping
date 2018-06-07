#!/bin/bash

#PBS -l walltime=24:00:00,mem=24gb,nodes=1:ppn=16
#PBS -N GWAS_stability_resampling
#PBS -M neyha001@umn.edu
#PBS -m abe
#PBS -r n

# Change the working directory
cd /panfs/roc/groups/6/smithkp/neyha001/QTLMapping/S2MET_Mapping/Scripts/GWAS

#module load R/3.4.0
module load R/3.2.0_intel_mkl

# Resampling GWAS
Rscript stability_mapping_resampling.R
