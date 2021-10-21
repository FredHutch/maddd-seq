#!/bin/bash

set -Eeuo pipefail

# Count up the number of aligned reads
samtools view -b -h --regions-file "${target_regions_bed}" "unmasked.bam" > "aligned.bam"
