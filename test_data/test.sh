#!/bin/bash

set -Eeuo pipefail

nextflow \
    run \
    -c ../nextflow.config \
    ../main.nf \
    -profile docker \
    --sample_sheet manifest.csv \
    --output output \
    -resume
