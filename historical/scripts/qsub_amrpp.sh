#!/bin/bash
#PBS -N test_serotyping
#PBS -l walltime=00:10:00
#PBS -l procs=16
#PBS -l pmem=1g
#PBS -q batch
#PBS -j oe
#PBS -A rakus

eval "$(conda shell.bash hook)"
conda activate nextflow
module load singularity

cd /mnt/home/groups/nmrl/bact_analysis/amrplusplus_v2/
nextflow run main_AmrPlusPlus_v2.nf -profile singularity --output test_results