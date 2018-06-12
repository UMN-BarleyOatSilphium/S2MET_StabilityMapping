#!/bin/bash

#PBS -l walltime=04:00:00,mem=24gb,nodes=1:ppn=1
#PBS -N marker_by_env_pop_param
#PBS -M neyha001@umn.edu
#PBS -m abe
#PBS -r n

# Change the working directory
cd /panfs/roc/groups/6/smithkp/neyha001/QTLMapping/S2MET_Mapping/Scripts/Marker_Stability/

# module load R/3.4.0
module load R/3.2.0_intel_mkl

# Marker effect by environment computation
Rscript marker_effect_by_env.R

