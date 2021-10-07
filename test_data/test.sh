#!/bin/bash

set -Eeuo pipefail

# The user must specimen a manifest
MANIFEST=$1

# Use the pre-compiled yeast genome by default
GENOME=${GENOME:-genome/yeast/GCF_000146045.2_R64_genomic.fna}

if (( ${#MANIFEST} == 0 )); then

    echo "Please specify manifest"

else

    echo "Running data in manifest: $MANIFEST"

    nextflow \
        run \
        -c ../nextflow.config \
        ../main.nf \
        -profile docker \
        --sample_sheet $MANIFEST \
        --output ${MANIFEST%.csv}.output \
        --genome ${GENOME}.'*' \
        -with-report ${MANIFEST%.csv}.output.html \
        -resume

fi