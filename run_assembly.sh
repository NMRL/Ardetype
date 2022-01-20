#!/bin/bash
#PBS -N test_assembly
#PBS -l walltime=24:00:00
#PBS -l procs=12
#PBS -l pmem=2g
#PBS -q batch
#PBS -j oe
#PBS -A rakus

module load singularity
cd ~/nmrl/bact_analysis/NMRL_Bact_Assembly_Inhouse/
singularity run snakemake.sif snakemake --keep-going --cores 12
