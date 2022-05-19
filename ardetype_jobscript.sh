#!/bin/bash
#PBS -N ardetype
#PBS -l walltime=24:00:00
#PBS -l procs=49
#PBS -l mem=512g
#PBS -q batch
#PBS -j oe
#PBS -A rakus

eval "$(conda shell.bash hook)"
DEFAULT_ENV=/mnt/home/$(whoami)/.conda/envs/mamba_env/envs/snakemake
conda activate $DEFAULT_ENV
snakefile=${1}
config_file=${2}

snakemake --snakefile ${snakefile} --configfile ${config_file} --keep-going --use-envmodules --use-conda --conda-frontend conda --rerun-incomplete --cores 12 --latency-wait 30 -np
#-np - dry run (for testing purposes)
#--keep-going - continue excecution if job fails
#--use-envmodules - to load singularity
#--conda-frontend conda - to use conda package manager (snakemake default is mamba)
#--rerun-incomplete - to avoid breaking while testing
#--forceall --rulegraph | dot -Tpdf > dag.pdf - to visualize jobs as graph in pdf format
#Attempts to run in cluster mode on snakemake 7.6.1 failed - need to run snakemake in cluster mode on login node 
