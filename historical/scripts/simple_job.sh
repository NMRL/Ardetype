#!/bin/bash
#PBS -N simple_job
#PBS -q batch
#PBS -l walltime=00:30:00
#PBS -l nodes=1:ppn=1
#PBS -j oe

echo "Hello from node `/bin/hostname`"

