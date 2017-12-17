#!/bin/bash

#PBS -l walltime=24:00:00,mem=24gb,nodes=1:ppn=16
#PBS -N S2_MET_qxe_mapping
#PBS -M neyha001@umn.edu
#PBS -m abe
#PBS -r n

# Change the working directory
cd /panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET_Mapping/Mapping

#module load R/3.4.0
module load R/3.2.0_intel_mkl

# Mapping of main effects and GxE
Rscript S2MET_mapping.R

# Mapping of stability
#Rscript S2MET_stability_mapping.R
