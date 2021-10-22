#!/bin/bash

set -Eeuo pipefail

# Use the pre-compiled yeast genome by default
GENOME=${GENOME:-genome/yeast/GCF_000146045.2_R64_genomic.fna.gz}

echo "Running on $GENOME"
nextflow \
    run \
    -c ../nextflow.config \
    ../index.nf \
    -profile docker \
    --output ${GENOME}/ \
    --genome_fasta ${GENOME} \
    -with-report index_test.output.html \
    -resume
