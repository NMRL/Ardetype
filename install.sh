#!/bin/bash
conda install -y -n base -c conda-forge mamba
mamba create -c conda-forge -c bioconda -n snakemake snakemake
conda env update -n snakemake --file environment.yaml