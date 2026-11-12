#!/bin/bash
#PBS -N ardetype
#PBS -l walltime=24:00:00
#PBS -l procs=49
#PBS -l mem=512g
#PBS -q batch
#PBS -j oe
#PBS -A rakus

eval "$(conda shell.bash hook)"
DEFAULT_ENV=$(dirname $(dirname $(which conda)))/envs/ardetype
conda activate $DEFAULT_ENV
snakefile=${1}
config_file=${2}

# Singularity selector (prefers system singularity, falls back to conda) must come before the conda env on PATH
export PATH="$(dirname $(dirname $(readlink -f ${snakefile})))/subscripts/bin:${PATH}"

snakemake --snakefile ${snakefile} --configfile ${config_file} --keep-going --use-envmodules --use-conda --conda-frontend conda --rerun-incomplete --cores 12 --latency-wait 30