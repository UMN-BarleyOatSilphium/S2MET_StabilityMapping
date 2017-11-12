#!/bin/bash

#PBS -l walltime=24:00:00,mem=24gb,nodes=1:ppn=16
#PBS -N S2_MET_gwas_simulation
#PBS -M neyha001@umn.edu
#PBS -m abe
#PBS -r n

# Change the working directory
cd /panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET_Mapping/Simulation

module load R/3.4.0

# Clustering and prediction
Rscript S2MET_mapping_simulation.R
