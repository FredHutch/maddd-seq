#!/bin/bash

set -Eeuo pipefail

echo "Specimen: $specimen"
echo "Reference genome index:"
echo "${ref}" | tr ' ' '\n'

# The input reads are organized into folders
for FOLDER in SSC_FWD_R1 SSC_FWD_R2 SSC_REV_R1 SSC_REV_R2; do

    echo \$FOLDER:
    ls -lahtr \$FOLDER
    echo

done

# The ref variable contains all of the index files
# To get the base filename, we will find the file ending
# with .amb and strip off that suffix
GENOME=\$(echo "${ref}" | tr ' ' '\\n' | grep '.amb' | sed 's/.amb//' )
echo "Reference genome prefix: \$GENOME"

# Run the alignment twice, once for each strand
echo """POS SSC_FWD_R1 SSC_REV_R2
NEG SSC_FWD_R2 SSC_REV_R1""" | while read STRAND FWD REV; do

    echo "Running BWA MEM for \$STRAND strand with \$FWD and \$REV"
    bwa \
        mem \
        -a \
        -t ${task.cpus} \
        -T ${params.min_align_score} \
        -C \
        "\$GENOME" \
        <(gunzip -c \$FWD/*) \
        <(gunzip -c \$REV/*) \
    | samtools \
        sort \
        --threads ${task.cpus} \
        -o \$STRAND.SSC.bam -

done

echo "DONE"