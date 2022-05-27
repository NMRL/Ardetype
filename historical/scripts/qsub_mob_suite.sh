#!/bin/bash
#PBS -N mob_suite
#PBS -l walltime=00:10:00
#PBS -l procs=12
#PBS -l pmem=1g
#PBS -q batch
#PBS -j oe
#PBS -A rakus

module load singularity
mob_sif=$(find /mnt/home/groups/nmrl/image_files/ -type f -name "mob_suite_3.0.3.sif")

contigs=${1}
sample_id=$(basename ${contigs})

singularity run ${mob_sif} mob_recon --infile ${contigs} --outdir ~/
singularity run ${mob_sif} mob_typer --csv ${contigs} >> ${output_path}${sample_id::-6}_legsta.csv