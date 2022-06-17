#!/bin/bash
reference_path=${1}
read_1=${2}
read_2=${3}
sample_id=$(echo $read_1 | grep -oP '[0-9]_S[0-9]_L001')
echo ${sample_id}
bwa index ${reference_path}
bwa mem ${reference_path} -v 3 ${read_1} ${read_2} | samtools view -bS - | samtools sort -o ${sample_id}_sorted.bam
samtools index ${sample_id}_sorted.bam
samtools coverage -A -w 100 ${sample_id}_sorted.bam > ${sample_id}_coverage.txt
# samtools depth -m 0 ${sample_id}_sorted.bam > ${sample_id}_seq_depth.txt
# python depth_plot.py ${sample_id}_seq_depth.txt ${sample_id}_seq_depth_plot.html