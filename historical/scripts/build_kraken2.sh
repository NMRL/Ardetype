#!/bin/bash
#PBS -N build_db
#PBS -l nodes=1:ppn=32,pmem=6g
#PBS -l walltime=24:00:00
#PBS -q long
#PBS -j oe
#PBS -A rakus

eval "$(conda shell.bash hook)"
conda activate kraken2

# cd ~/enterovirus/vipr_refseq_enteroviruses/
# kraken2-build --add-to-library 62231528082-GenomicFastaResults_kraken2.fasta --db vipr_enterovir
# kraken2-build --build --db vipr_enterovir --threads 32
kraken2-build --download-taxonomy --db viral --use-ftp
kraken2-build --download-library viral --db viral --use-ftp
kraken2-build --build --db viral --threads 32