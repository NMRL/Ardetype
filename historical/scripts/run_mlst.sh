#!/bin/bash -i
input=${1}
output=${2}
image=${3}

module load singularity
singularity run ${image} mlst --csv ${input} >> ${output}