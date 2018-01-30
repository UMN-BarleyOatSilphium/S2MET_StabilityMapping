#!/bin/bash

#PBS -l walltime=24:00:00,mem=62gb,nodes=1:ppn=24
#PBS -N S2_MET_stability_mxe
#PBS -M neyha001@umn.edu
#PBS -m abe
#PBS -r n

# Change the working directory
cd /panfs/roc/groups/6/smithkp/neyha001/Genomic_Selection/S2MET_Mapping/PhenotypeManipulation

# module load R/3.4.0
module load R/3.2.0_intel_mkl

# Marker effect by environment computation
Rscript S2MET_marker_effect_by_env.R HeadingDate
Rscript S2MET_marker_effect_by_env.R GrainYield
Rscript S2MET_marker_effect_by_env.R PlantHeight