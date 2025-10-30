#!/bin/bash

set -e

echo
echo "Running trim_overhang.py"
echo

ls -lahtr
echo

# Use a random string to keep the files paired
RAND_ID="\$(cat /dev/urandom | tr -dc "A-Za-z" | head -c 8)"

trim_overhang.py \
    --input-bam "untrimmed.*.bam" \
    --input-positions read_positions.csv.gz \
    --output-read1 "${specimen}_${shard_ix}_\${RAND_ID}_R1.fastq.gz" \
    --output-read2 "${specimen}_${shard_ix}_\${RAND_ID}_R2.fastq.gz"

echo
echo DONE
echo



ls -lahtr