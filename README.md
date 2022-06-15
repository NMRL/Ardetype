# ARDETYPE (v.1.0)

NGS data processing pipeline designed to perform species-agnostic and species-specific analysis of short paired-end (PE) bacterial reads.

Pipeline is structured in terms of modules. Each module corresponds to a Snakemake script (snakefile) and a python module object. Snakefiles are used to define data processing rules. 

Module class objects store information about inputs and expected outputs for these rules, handle file placement operations, module configuration, information transfer between modules, and other processes that are not directly involved in running analysis on files.

|Module      |Type     |Description                                                         |Tools                                             |
|:----------:|--------:|:-------------------------------------------------------------------|:------------------------------------------------:|
|bact_core   |agnostic |QC, host filtering, denovo assembly, taxonomic classification       |[fastp](https://github.com/OpenGene/fastp), [kraken2](https://github.com/DerrickWood/kraken2), [shovill](https://github.com/tseemann/shovill), [krakentools](https://github.com/jenniferlu717/KrakenTools), [krona](https://github.com/marbl/Krona)       |
|bact_shell  |agnostic |assembly QC, resistance profiling, plasmid reconstruction & typing  |[quast](https://github.com/ablab/quast), [rgi-card](https://github.com/arpcard/rgi), [amr++v2.0](https://megares.meglab.org/amrplusplus/latest/html/v2/), [resfinder](https://bitbucket.org/genomicepidemiology/resfinder/src/master/), [mob-suite](https://github.com/phac-nml/mob-suite)  |
|bact_tip    |specific |species-dependent sub-typing                                        |[hicap](https://github.com/scwatts/hicap), [meningotype](https://github.com/MDU-PHL/meningotype), [legsta](https://github.com/tseemann/legsta), [Kleborate](https://github.com/katholt/Kleborate/wiki), [AgrVATE](https://github.com/VishnuRaghuram94/AgrVATE), [spaTyper]( https://github.com/HCGB-IGTP/spaTyper), [Staphopia-sccmec](https://github.com/staphopia/staphopia-sccmec), [emmtyper](https://github.com/MDU-PHL/emmtyper), [seqsero](https://github.com/denglab/SeqSero), [sistr](https://github.com/phac-nml/sistr_cmd), [lissero](https://github.com/MDU-PHL/LisSero), [PubMLST database API](https://pubmlst.org/), [Institute Pasteur MLST database API](https://bigsdb.pasteur.fr/), [Legionella pneumophila in silico Serogroup Prediction](https://github.com/NMRL/legionella_pneumophila_genomics), [ectyper](https://github.com/phac-nml/ecoli_serotyping), [seroba](https://github.com/sanger-pathogens/seroba)|

## Configuration
Pipeline is designed to be run by NMRL users on RTU HPC, where HPC-level configuration is available out-of-the-box. 
|Configuration level|Dependencies|
|-------------------|------------|
|HPC                |Torque/PBS, Conda, Singularity          |
|Conda              |Snakemake (should be installed for each user using -s flag to the ardetype.py script)            |
|Python             |numpy==1.22.3, pandas==1.4.2, PyYAML==6.0, requests==2.27.1, bs4==0.0.1           |
|Kraken2            |[Pre-built](https://ccb.jhu.edu/software/kraken2/downloads.shtml) or custom databases for human and bacteria|
|Resfinder4         |[Database](https://bitbucket.org/genomicepidemiology/resfinder_db/src)|

## Installation
To install from scratch, you will need an access to a Linux system with **root access** and installed [singularity](https://sylabs.io/guides/3.0/user-guide/installation.html) to build containers. [WSL](https://docs.microsoft.com/en-us/windows/wsl/install) or [Virtual Machine](https://www.arcserve.com/blog/dead-simple-guide-installing-linux-virtual-machine-windows) should also work. 

Clone the repository to your local machine and use [singularity recipe files](https://github.com/NMRL/NMRL_Bact_Assembly_Inhouse/tree/ardetype/config_files/s_recipes) to build containers,<br>then copy to HPC cluster so that they can be accessed by the pipeline scripts.

Clone the repository to the cluster and edit files found in config_files folder to match your local setup:
|File|Scope|
|----|-----|
|module_data.json|paths to cluster_config file and snakefiles|
|config_modular.yaml|paths to singularity image files, [kraken2](https://github.com/DerrickWood/kraken2) databases, [resfinder](https://bitbucket.org/genomicepidemiology/resfinder/src/master/) [database](https://bitbucket.org/genomicepidemiology/resfinder_db/src), path to [Legionella pneumophila in silico Serogroup Prediction](https://github.com/NMRL/legionella_pneumophila_genomics) tool|

## Using the pipeline
 - Note: pipeline accepts only fastq files that are named according to illumina conventions (sample_id_R{1,2}_001.fastq.gz).
 - Testing (to see what jobs will be executed): 
     
     ``` python ardetype.py -t -i path_to_folder_with_fastq/ -o path_to_output_folder -m all ```
 
 - Running: 
     
     ``` python ardetype.py -i path_to_folder_with_fastq/ -o path_to_output_folder -m all ```
 
