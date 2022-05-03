#!/bin/bash
#PBS -N ardetype
#PBS -l walltime=24:00:00
#PBS -l procs=48
#PBS -l mem=256g
#PBS -q batch
#PBS -j oe
#PBS -A rakus

eval "$(conda shell.bash hook)"
conda activate /mnt/home/jevgen01/.conda/envs/mamba_env/envs/snakemake
snakefile = ${1}
config_file = ${2}
snakemake --snakefile ${snakefile} --keep-going --use-envmodules --use-conda --conda-frontend conda --rerun-incomplete --cores 48 --latency-wait 30
#--keep-going - continue excecution if job fails
#--use-envmodules - to load singularity
#--conda-frontend conda - to use conda package manager (snakemake default is mamba)
#--rerun-incomplete - to avoid breaking while testing