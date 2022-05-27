#!/bin/bash
#PBS -N test_assembly
#PBS -l walltime=01:00:00
#PBS -l procs=16
#PBS -l pmem=1g
#PBS -q batch
#PBS -j oe
#PBS -A rakus

eval "$(conda shell.bash hook)"
conda activate rgi
module load singularity

fastq_sif_path=$(find /mnt/home/groups/nmrl/image_files/ -type f -name "fastq_processing.sif")
qualimap_sif_path=$(find /mnt/home/groups/nmrl/image_files/ -type f -name "qualimap_latest.sif")
mlst_sif_path=$(find /mnt/home/groups/nmrl/bact_analysis/NMRL_Bact_Assembly_Inhouse/ -type f -name "mlst_quast.sif")
rgi_prokkka_sif_path=$(find /mnt/home/groups/nmrl/bact_analysis/NMRL_Bact_Assembly_Inhouse/ -type f -name "rgi_prokka.sif")

contigs=${1}

#Prokka
singularity run ${rgi_prokkka_sif_path} prokka --outdir ~/${contigs::-6}_prokka/ ${contigs}

#rgi
rgi main --input_sequence ${contigs} --output_file ~/${contigs::-6}.rgi --input_type contig --clean

#mlst
singularity run ${mlst_sif_path} mlst --csv ${contigs} >> ~/${contigs::-6}_mlst.csv