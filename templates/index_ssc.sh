#!/bin/bash

set -euo pipefail

# Get the output filepath from the channel
OUTPUT_BAM="${output_filename}"

# Sort the BAM
samtools sort -o "\$OUTPUT_BAM" unfiltered.bam
echo "Sorting unfiltered.bam -> \$OUTPUT_BAM"

# Indexing the BAM
echo "Indexing \$OUTPUT_BAM"
samtools index "\$OUTPUT_BAM"
