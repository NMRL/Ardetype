#!/bin/bash
#PBS -N resfinder
#PBS -l walltime=00:10:00
#PBS -l procs=12
#PBS -l pmem=1g
#PBS -q batch
#PBS -j oe
#PBS -A rakus

module load singularity

read_1=${1}
read_1_file=$(basename ${read_1})
read_2=${2}
read_2_file=$(basename ${read_2})
output=${3}
resfinder_image=/mnt/home/groups/nmrl/image_files/resfinder.sif

cp -r /mnt/home/groups/nmrl/db/db-resfinder/ ~/

cp ${read_1} ~/
cp ${read_2} ~/
singularity run ${resfinder_image} run_resfinder.py -ifq ~/${read_1_file} ~/${read_2_file} -acq -l 0.6 -t 0.9 -l_p 0.6 -t_p 0.9 -db_res ~/db-resfinder/db_resfinder -o ~/resfinder_output
mv ~/resfinder_output/* ${output}
rm -r ~/db-resfinder ~/${read_1_file} ~/${read_2_file} ~/resfinder_output


#https://cge.food.dtu.dk/services/ResFinder/instructions.php

#   -h, --help            show this help message and exit
#   -ifa INPUTFASTA, --inputfasta INPUTFASTA
#                         Input fasta file.
#   -ifq INPUTFASTQ [INPUTFASTQ ...], --inputfastq INPUTFASTQ [INPUTFASTQ ...]
#                         Input fastq file(s). Assumed to be single-end fastq if only one file is provided, and assumed to be paired-end data if two files are provided.
#   -o OUT_PATH, --outputPath OUT_PATH
#                         Path to blast output
#   -b BLAST_PATH, --blastPath BLAST_PATH
#                         Path to blastn
#   -k KMA_PATH, --kmaPath KMA_PATH
#                         Path to KMA
#   -s SPECIES, --species SPECIES
#                         Species in the sample
#   -db_res DB_PATH_RES, --db_path_res DB_PATH_RES
#                         Path to the databases for ResFinder
#   -db_res_kma DB_PATH_RES_KMA, --db_path_res_kma DB_PATH_RES_KMA
#                         Path to the ResFinder databases indexed with KMA. Defaults to the 'kma_indexing' directory inside the given database directory.
#   -d DATABASES, --databases DATABASES
#                         Databases chosen to search in - if none is specified all is used
#   -acq, --acquired      Run resfinder for acquired resistance genes
#   -ao ACQ_OVERLAP, --acq_overlap ACQ_OVERLAP
#                         Genes are allowed to overlap this number of nucleotides. Default: 30.
#   -l MIN_COV, --min_cov MIN_COV
#                         Minimum (breadth-of) coverage of ResFinder within the range 0-1.
#   -t THRESHOLD, --threshold THRESHOLD
#                         Threshold for identity of ResFinder within the range 0-1.
#   -nano, --nanopore     If nanopore data is used
#   -c, --point           Run pointfinder for chromosomal mutations
#   -db_point DB_PATH_POINT, --db_path_point DB_PATH_POINT
#                         Path to the databases for PointFinder
#   -db_point_kma DB_PATH_POINT_KMA, --db_path_point_kma DB_PATH_POINT_KMA
#                         Path to the PointFinder databases indexed with KMA. Defaults to the 'kma_indexing' directory inside the given database directory.
#   -g SPECIFIC_GENE [SPECIFIC_GENE ...]
#                         Specify genes existing in the database to search for - if none is specified all genes are included in the search.
#   -u, --unknown_mut     Show all mutations found even if in unknown to the resistance database
#   -l_p MIN_COV_POINT, --min_cov_point MIN_COV_POINT
#                         Minimum (breadth-of) coverage of Pointfinder within the range 0-1. If None is selected, the minimum coverage of ResFinder will be used.
#   -t_p THRESHOLD_POINT, --threshold_point THRESHOLD_POINT
#                         Threshold for identity of Pointfinder within the range 0-1. If None is selected, the minimum coverage of ResFinder will be used.
#   -v, --version         Show program's version number and exit
#   --pickle              Create a pickle dump of the Isolate object. Currently needed in the CGE webserver. Dependency and this option is being removed.