#!/bin/bash

reference_fasta=${1}
read_1=${2}
read_1_name=${read_1::-9}
read_2=${3}
read_2_name=${read_2::-9}

bwa index $reference_fasta
bwa mem $reference_fasta -v 3 $read_1 $read_2 | samtools view -bS - | samtools sort -o ${read_1_name}_reference_based.bam
bedtools bamtofastq -i ${read_1_name}_reference_based.bam -fq ${read_1_name}_reference_based.fastq -fq2 ${read_2_name}_reference_based.fastq
gzip ${read_1_name}_reference_based.bam
gzip ${read_2_name}_reference_based.fastq
