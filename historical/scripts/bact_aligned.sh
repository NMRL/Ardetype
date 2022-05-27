#!/bin/bash
#PBS -N test_alignment
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

read_1=${1}
read_2=${2}
reference=${3}

# align to reference with bwa to bam
# index reference
singularity run $fastq_sif_path bwa index ${reference}
# align reads
singularity run $fastq_sif_path bwa mem ${reference} -v 3 ${read_1} ${read_2}| singularity run $fastq_sif_path samtools view -bS - | singularity run $fastq_sif_path samtools sort -o ~/${read_1::-9}_sorted.bam

# produce vcf from bam
singularity run $fastq_sif_path freebayes -f ${reference} ~/${read_1::-9}_sorted.bam | singularity run $fastq_sif_path vcffilter -f "( QUAL > 20 )" AND "( DP > 25 )" > ~/${read_1::-9}.vcf

#generate consensus with bcftools 
singularity run $fastq_sif_path bgzip -c ~/${read_1::-9}.vcf > ~/${read_1::-9}.vcf.gz
singularity run $fastq_sif_path tabix -p vcf ~/${read_1::-9}.vcf.gz
cat ${reference} | singularity run $fastq_sif_path bcftools consensus ~/${read_1::-9}.vcf.gz > ~/${read_1::-9}_consensus.fa

#Prokka
singularity run ${rgi_prokkka_sif_path} prokka --outdir ~/${read_1::-9}_prokka/ ~/${read_1::-9}_consensus.fa

#rgi
rgi main --input_sequence ~/${read_1::-9}_consensus.fa --output_file ~/${read_1::-9}.rgi --input_type contig --clean

#mlst
singularity run ${mlst_sif_path} mlst --csv ~/${read_1::-9}_consensus.fa >> ~/${read_1::-9}_mlst.csv