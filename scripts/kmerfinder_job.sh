#!/bin/bash
#PBS -N kfinder
#PBS -l walltime=0:20:00
#PBS -l procs=12
#PBS -l pmem=2g
#PBS -q batch
#PBS -j oe
#PBS -A rakus

module load singularity
kmerfinder_sif_path="kmerfinder_3.0.2.sif"
file_1_path=${1}
file_2_path=${2}
sample_id=${3}
db=${4}

file_1=$(basename ${file_1_path})
file_2=$(basename ${file_2_path})


if [[ $db == 'viral' ]]; then
    db_path=db-kmerfinder/${db}/${db}.TG
else
    db_path=db-kmerfinder/${db}/${db}.ATG
fi

cd /mnt/home/groups/nmrl/bact_analysis/kmerfinder_suite
mkdir -p output/${sample_id}
singularity run ${kmerfinder_sif_path} -i input_files/${file_1} input_files/${file_2} -db ${db_path} -o output/${sample_id}