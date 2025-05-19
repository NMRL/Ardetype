#!/bin/bash

# Initialize conda shell environment
eval "$(conda shell.bash hook)"

# Arguments
TOOL_NAME=$1
TOOL_VERSION=$2

# Dynamic naming logic
TAG="${TOOL_NAME}_${TOOL_VERSION}"

# Configurations
CONDA_CHANNELS="-c bioconda -c conda-forge"
ENV_DEF_PATH="resources/env_defs/"         
RECIPE_PATH="resources/conda_recipes/"
IMAGE_PATH="resources/image_files/${TAG}.sif"

# Options
BUILD_CONDA_ENV=false
CHECKSUM_CONDA_DEF=false
CLEAR_CONDA_ENV=false
CLEAR_IMAGE=false
CLEAR_ENV_DEF=true
CLEAR_RECIPE=true
FETCH_ENV_DEF=true
CHECKSUM_IMAGE=false
BUILD_SIF=true

# Cleanup existing Image if needed
if $CLEAR_IMAGE && [ -f "$IMAGE_PATH" ]; then
    rm "$IMAGE_PATH"
    rm "${IMAGE_PATH}.sha"    
else 
    echo "Image file for $TAG will not be removed."
fi

# Cleanup Conda Environment if needed
if $CLEAR_CONDA_ENV; then
    if conda env list | grep -q "$TAG"; then
        echo "Environment $TAG exists. Removing it..."
        conda env remove -n "$TAG" --yes
        rm "${ENV_DEF_PATH}${TAG}.yaml.sha"
    else 
        echo "Environment $TAG does not exist."
    fi
    else 
        echo "Conda environment for $TAG will not be removed."
fi

# If environment creation is needed
if $BUILD_CONDA_ENV; then
    conda create -n "$TAG" --yes  # Create the environment after removal
    conda activate "$TAG"
    mamba install ${CONDA_CHANNELS} ${TOOL_NAME}=${TOOL_VERSION} --yes
    conda env export > "${ENV_DEF_PATH}${TAG}.yaml"
    # Clean YAML (remove last line)
    sed -i '$d' "${ENV_DEF_PATH}${TAG}.yaml"
    conda deactivate

    # Checksum conda environment file
    if $CHECKSUM_CONDA_DEF && [ -f "${ENV_DEF_PATH}${TAG}.yaml" ]; then
        sha256sum "${ENV_DEF_PATH}${TAG}.yaml" > "${ENV_DEF_PATH}${TAG}.yaml.sha"
    else 
        echo "Checksum on environment definition file for $TAG will not be generated."
    fi
elif $FETCH_ENV_DEF; then
    # Get the environment definition
    python get_conda_def.py --package $TOOL_NAME --version $TOOL_VERSION --output "${ENV_DEF_PATH}${TAG}.yaml"
    # Checksum conda environment file
    if $CHECKSUM_CONDA_DEF && [ -f "${ENV_DEF_PATH}${TAG}.yaml" ]; then
        sha256sum "${ENV_DEF_PATH}${TAG}.yaml" > "${ENV_DEF_PATH}${TAG}.yaml.sha"
    else 
        echo "Checksum on environment definition file for $TAG will not be generated."
    fi
else
    echo "Environment definition file for $TAG will not be created."
fi

# Create Singularity Recipe
cat <<EOF > "${RECIPE_PATH}${TAG}.recipe"
Bootstrap: docker
From: continuumio/miniconda3

%files
    ../env_defs/${TAG}.yaml ${TAG}.yaml
%post
    conda env create --name ${TAG} --file ${TAG}.yaml
%environment
    export PATH=/opt/conda/envs/${TAG}/bin:\$PATH
%runscript
    exec "\$@"
EOF

# Build the Singularity Image
if $BUILD_SIF ; then
    cd "$(dirname $IMAGE_PATH)"
    singularity build --fakeroot "${TAG}.sif" "../conda_recipes/${TAG}.recipe"
else
    echo "Singularity image for $TAG will not be created."
fi

# Cleanup Environment Definition if needed
if $CLEAR_ENV_DEF && [ -f "${ENV_DEF_PATH}${TAG}.yaml" ]; then
    rm "${ENV_DEF_PATH}${TAG}.yaml"
    rm "${ENV_DEF_PATH}${TAG}.yaml.sha"
else
    echo "Environment definition file for $TAG will not be removed."
fi

# Cleanup Recipe if needed
if $CLEAR_RECIPE && [ -f "${RECIPE_PATH}${TAG}.recipe" ]; then
    rm "${RECIPE_PATH}${TAG}.recipe"
else
    echo "Recipe file for $TAG will not be removed."
fi

if $CHECKSUM_IMAGE; then
    cd ..
    CHECKSUM_PATH=${IMAGE_PATH}.sha
    sha256sum ${IMAGE_PATH} > $CHECKSUM_PATH
else
    echo "sha256 checksum on image file for $TAG will not be generated."
fi
