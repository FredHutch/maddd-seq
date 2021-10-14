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
    -o "${R1.name.replaceAll(/.fastq.gz/, '')}.trimmed.fastq.gz" \
    -p "${R2.name.replaceAll(/.fastq.gz/, '')}.trimmed.fastq.gz" \
    --json="${specimen}.cutadapt.json" \
    --cores=${task.cpus} \
    -a "${params.RD1_ADAPTER_3P}" \
    -A "${params.RD2_ADAPTER_3P}" \
    "$R1" \
    "$R2"

echo DONE
