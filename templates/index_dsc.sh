#!/bin/bash

set -Eeuo pipefail

# Sort the DSC alignments
samtools sort -o DSC.bam unsorted.DSC.bam

# Make an index file
samtools index DSC.bam