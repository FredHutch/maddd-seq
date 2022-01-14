#!/bin/bash

set -e

echo
echo "Running trim_overhang.py"
echo

ls -lahtr
echo

trim_overhang.py \
    --input-bam untrimmed.bam \
    --input-positions read_positions.csv.gz \
    --output-read1 "${specimen}_${shard_ix}_R1.fastq.gz" \
    --output-read2 "${specimen}_${shard_ix}_R2.fastq.gz"

echo
echo DONE
echo



ls -lahtr