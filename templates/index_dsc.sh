#!/bin/bash

set -Eeuo pipefail

# Sort the DSC alignments
samtools sort --threads ${task.cpus} -o DSC.bam unsorted.DSC.bam

# Make an index file
samtools index DSC.bam