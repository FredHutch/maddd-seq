#!/bin/bash

set -Eeuo pipefail

nextflow \
    run \
    ../main.nf \
    --manifest manifest.csv \
    --output output
