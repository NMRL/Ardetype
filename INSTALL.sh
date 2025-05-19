#!/bin/bash

# Create images from anaconda packages
CONDA_LIST="resources/anaconda_list.csv"
SCRIPT_PATH="build_image.sh"
PARALLEL_JOBS=4

export SCRIPT_PATH  # So it's visible in the subshells

cat "$CONDA_LIST" | parallel -j $PARALLEL_JOBS --colsep ',' '
    echo "Running: $SCRIPT_PATH {1} {2}"
    bash "$SCRIPT_PATH" "{1}" "{2}"
'

# Create images from custom recipes
CUSTOM_LIST="resources/custom_list.csv"

cat "$CUSTOM_LIST" | parallel -j $PARALLEL_JOBS --colsep ',' '
    echo "Building from custom image: {1}";
    singularity build --fakeroot {2} {1};
'

# Create images from dockerhub
CUSTOM_LIST="resources/dockerhub_list.csv"

cat "$CUSTOM_LIST" | parallel -j $PARALLEL_JOBS --colsep ',' '
    echo "Pulling official image from hub.docker.com: {1} {2}";
    cd resources/image_files ;
    singularity pull docker://{2} ;
    cd ../../;
'

# # Download databases
# 1. T>T human
# 2. phix size change
# 3. Compress
# /mnt/beegfs2/home/groups/nmrl/db/T2T-CHM13v2
# Initialize config

