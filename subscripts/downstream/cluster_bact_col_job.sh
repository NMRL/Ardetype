#!/bin/bash
#PBS -N clst_bct_cl
#PBS -l walltime=03:30:00
#PBS -l procs=32
#PBS -l mem=64g
#PBS -q batch
#PBS -j oe
#PBS -A rakus

SCRIPT_PATH=/mnt/beegfs2/home/groups/nmrl/bact_analysis/Ardetype/subscripts/downstream

cd $SCRIPT_PATH && python cluster_bact_collections.py