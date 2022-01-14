#!/bin/bash

set -euo pipefail

# Copy the header
samtools view -H ${bam} > adduct.reads.sam

# Add the reads from the families which yielded adducts
samtools view ${bam} \
    | fgrep -f <(fgrep -f <(gunzip -c ${adduct_families}) <(gunzip -c ${read_families}) | sed 's/,.*//') \
    >> adduct.reads.sam

# Convert SAM -> BAM
samtools view -bS adduct.reads.sam > adduct.reads.bam

rm adduct.reads.sam