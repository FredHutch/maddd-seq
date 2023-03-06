#!/bin/bash

set -euo

# Copy the header
echo "Copying header"
samtools view -H ${bam} > adduct.reads.sam

# Add the reads from the families which yielded adducts
echo "Getting families from this shard which contain adducts"
( fgrep -f <(gunzip -c ${adduct_families}) <(gunzip -c ${read_families}) > shard.adduct.families.txt  || test $? = 1 )

# If there are families from this shard which do contain adducts
if (( \$(cat shard.adduct.families.txt | wc -l) > 0 )); then
    echo "Extracting reads"

    # Note that this step _will_ raise an error if no matches are found
    samtools view ${bam} \
        | fgrep -f <(cat shard.adduct.families.txt | sed 's/,.*//') \
        >> adduct.reads.sam
fi

# Convert SAM -> BAM
echo "Converting SAM to BAM"
samtools view -bS adduct.reads.sam > adduct.reads.bam

echo "Removing temporary SAM file"
rm adduct.reads.sam

echo "Done"