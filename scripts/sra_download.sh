#!/bin/bash
#PBS -N download_reads
#PBS -l nodes=1:ppn=24,pmem=2g
#PBS -l walltime=24:00:00
#PBS -q long
#PBS -j oe
#PBS -A rakus

eval "$(conda shell.bash hook)"
conda activate bio_env

cd ~/sra
cat ~/systemic_salmonella_list.csv | xargs -n1 -P23 -I% fastq-dump --split-3 --gzip --accession %
    