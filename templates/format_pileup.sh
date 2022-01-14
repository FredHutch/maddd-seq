#!/bin/bash

set -Eeuo pipefail

echo "Specimen: $specimen"
echo "Reference genome index:"
echo "${ref}" | tr ' ' '\n'

# The ref variable contains all of the index files
# To get the base filename, we will find the file ending
# with .amb and strip off that suffix
GENOME=\$(echo "${ref}" | tr ' ' '\\n' | grep '.amb' | sed 's/.amb//' )
echo "Reference genome prefix: \$GENOME"

# If there is no uncompressed genome file
if [[ ! -s \$GENOME ]]; then

    # There must be a gzip-compressed file
    [[ ! -s \${GENOME}.gz ]] && "Cannot find gzip-compressed genome reference"
    # which we can use to make the uncompressed version
    gunzip -c \${GENOME}.gz > \${GENOME}

fi

# Make a pileup file from the DSC data
samtools \
    mpileup \
    -f \${GENOME} \
    DSC.bam \
| gzip -c \
> DSC.pileup.gz