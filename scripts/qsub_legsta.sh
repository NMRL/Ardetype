#!/bin/bash
#PBS -N legsta
#PBS -l walltime=00:10:00
#PBS -l procs=12
#PBS -l pmem=1g
#PBS -q batch
#PBS -j oe
#PBS -A rakus

module load singularity

legsta_sif=$(find /mnt/home/groups/nmrl/image_files/ -type f -name "legsta_latest.sif")

contigs=${1}
output_path=${2}
sample_id=$(basename ${contigs})

singularity run ${legsta_sif} legsta --csv ${contigs} >> ${output_path}${sample_id::-6}_legsta.csv