#!/bin/bash
#PBS -N test_assembly
#PBS -l walltime=24:00:00
#PBS -l procs=12
#PBS -l pmem=2g
#PBS -q batch
#PBS -j oe
#PBS -A rakus

eval "$(conda shell.bash hook)" 
conda activate /mnt/home/jevgen01/.conda/envs/mamba_env/envs/snakemake
cd ~/nmrl/bact_analysis/NMRL_Bact_Assembly_Inhouse/
snakemake --keep-going --use-envmodules --cores 12
