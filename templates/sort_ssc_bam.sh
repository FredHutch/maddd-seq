#!/bin/bash

set -Eeuo pipefail

# Sort and index both of the SSC BAM files
for INPUT_FP in unsorted.*; do

    # Remove the prefix
    OUTPUT_FP=\${INPUT_FP#unsorted.}

    # Sort the BAM
    samtools sort --threads ${task.cpus} -o \$OUTPUT_FP \$INPUT_FP

    # Make an index file
    samtools index \$OUTPUT_FP

done
