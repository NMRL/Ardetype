#!/bin/bash
#PBS -N test_serotyping
#PBS -l walltime=00:10:00
#PBS -l procs=16
#PBS -l pmem=1g
#PBS -q batch
#PBS -j oe
#PBS -A rakus

eval "$(conda shell.bash hook)"
# conda activate rgi
module load singularity

mlst_sif_path=$(find /mnt/home/groups/nmrl/bact_analysis/NMRL_Bact_Assembly_Inhouse/ -type f -name "mlst_quast.sif")
rgi_prokkka_sif_path=$(find /mnt/home/groups/nmrl/bact_analysis/NMRL_Bact_Assembly_Inhouse/ -type f -name "rgi_prokka.sif")
sistr_path=$(find /mnt/home/groups/nmrl/image_files/ -type f -name "sistr_latest.sif")
seqsero2_path=$(find /mnt/home/groups/nmrl/image_files/ -type f -name "seqsero2_latest.sif")
hamr_path=$(find /mnt/home/groups/nmrl/image_files/ -type f -name "hamronization_latest.sif")

contigs=${1}
read_1=${2}
read_2=${3}
sample_id=$(basename ${contigs})

mv ${read_1} ~/test/
mv ${read_2} ~/test/
mv ${contigs} ~/test/
# #rgi
# rgi main --input_sequence ${contigs} --output_file ~/salmonella_sra/rgi/${sample_id::-6}.rgi --input_type contig --clean

# #hamronize rgi
# mkdir ~/salmonella_sra/rgi/hamronized/
# for i in ~/salmonella_sra/rgi/*.txt; do singularity run ${hamr_path} rgi --input_file_name $(basename $i) --analysis_software_version rgi_v5 --analysis_software_version rgi_v1 --reference_database_version card_v1 ${i} --format json > ~/salmonella_sra/rgi/hamronized/$(basename ${i})_harm.json ; done
# singularity run ${hamr_path} summarize -o ~/salmonella_sra/rgi_summary.html -t interactive ~/salmonella_sra/rgi/hamronized/*.json

# #mlst
# singularity run ${mlst_sif_path} mlst --csv ${contigs} >> ~/salmonella_sra/mlst/${sample_id::-6}_mlst.csv

#sistr
singularity run ${sistr_path} sistr --qc -f csv -o ~/test/sistr/${sample_id::-6}_sistr.csv ~/test/${sample_id}
# --qc - check data quality - informative messages regarding missing loci, genes, issues with genome sizes or else
# -f output format

#seqsero2
singularity run ${seqsero2_path} SeqSero2_package.py -d ~/test/seqsero2/${sample_id::-6} -n ${sample_id::-6} -p 16 -t 2 -i ~/test/$(basename ${read_1}) ~/test/$(basename ${read_2})
# -t paired-end reads
# -n sample name
# -p threads
# -d directory name

mv ~/test/$(basename ${read_1}) ${read_1}
mv ~/test/$(basename ${read_2}) ${read_2}
mv ~/test/$(basename ${contigs}) ${contigs}