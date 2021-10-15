#!/bin/bash

set -Eeuo pipefail

echo "Specimen: $specimen"
echo "R1: $R1"
echo "R2: $R2"
echo "Reference genome index:"
echo "${ref}" | tr ' ' '\n'

# The ref variable contains all of the index files
# To get the base filename, we will find the file ending
# with .amb and strip off that suffix
GENOME=\$(echo "${ref}" | tr ' ' '\\n' | grep '.amb' | sed 's/.amb//' )
echo "Reference genome prefix: \$GENOME"

echo "Running BWA MEM, filtering out unmapped and secondary alignments, and sorting"
bwa \
    mem \
    -t ${task.cpus} \
    -T ${params.min_align_score} \
    -C \
    "\$GENOME" \
    ${R1} \
    ${R2} | \
samtools \
    view \
    -F 722 \
    -h \
    -b | \
samtools \
    sort \
    --threads ${task.cpus} \
    -o \
    aligned.bam

echo "DONE"