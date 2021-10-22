#!/bin/bash

set -Eeuo pipefail

# Count up the number of aligned reads
samtools view --threads ${task.cpus} -b -h --regions-file "${target_regions_bed}" "unmasked.bam" > "aligned.bam"
