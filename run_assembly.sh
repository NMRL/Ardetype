#!/bin/bash
#PBS -N test_assembly
#PBS -l walltime=01:00:00
#PBS -l procs=12
#PBS -l mem=24g
#PBS -q batch
#PBS -j oe
#PBS -A rakus

eval "$(conda shell.bash hook)"
conda activate /mnt/home/jevgen01/.conda/envs/mamba_env/envs/snakemake
cd ~/nmrl/bact_analysis/NMRL_Bact_Assembly_Inhouse/
snakemake --keep-going --use-envmodules --use-conda --conda-frontend conda --rerun-incomplete --cores 48 --latency-wait 30
#--keep-going - continue excecution if job fails
#--use-envmodules - to load singularity
#--conda-frontend conda - to use conda package manager (snakemake default is mamba)
#--rerun-incomplete - to avoid breaking while testing