#!/bin/bash

set -Eeuo pipefail

# The user must specimen a manifest
MANIFEST=$1

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
        -with-report ${MANIFEST%.csv}.output.html \
        -resume

fi