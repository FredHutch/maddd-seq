#!/bin/bash

set -Eeuo pipefail

echo "Processing specimen: $specimen"
echo "R1: $R1"
echo "R2: $R2"
echo "--quality-cutoff=${params.min_qvalue}"
echo "--minimum-length=${params.min_align_score}"
echo "--json="${specimen}.cutadapt.json""

cutadapt \
    --pair-filter=any \
    --quality-cutoff=${params.min_qvalue} \
    --minimum-length=${params.min_align_score} \
    -o trimmed.R1.fastq.gz \
    -p trimmed.R1.fastq.gz \
    --json="${specimen}.cutadapt.json" \
    "$R1" \
    "$R2"

echo DONE
