#!/bin/bash

set -Eeuo pipefail

echo "Processing specimen: $specimen"
echo "R1: $R1"
echo "R2: $R2"
echo "--cut=${params.barcode_length}"
echo "--minimum-length=${params.min_align_score}"
echo "--json=${specimen}.cutadapt.json"
echo "--rename='{id} {r1.cut_prefix}{r2.cut_prefix}'"

cutadapt \
    --pair-filter=any \
    --cut=${params.barcode_length} \
    -U ${params.barcode_length} \
    --minimum-length=${params.min_align_score} \
    -o "${R1.name.replaceAll(/.fastq.gz/, '')}.clipped.fastq.gz" \
    -p "${R2.name.replaceAll(/.fastq.gz/, '')}.clipped.fastq.gz" \
    --json="${specimen}.cutadapt.json" \
    --rename='{id} BC:Z:{r1.cut_prefix}{r2.cut_prefix}' \
    --cores ${task.cpus} \
    "$R1" \
    "$R2"

# Count up the number of times that each barcode was seen
echo "count,barcode" > barcode_counts.csv
gunzip -c "${R1.name.replaceAll(/.fastq.gz/, '')}.clipped.fastq.gz" \
    | awk 'NR % 4 == 1' \
    | sed 's/.* //' \
    | sort \
    | uniq -c \
    | sort -nrk1 \
    | sed 's/^ *//' \
    | tr ' ' ',' \
    >> barcode_counts.csv

gzip barcode_counts.csv

echo DONE
