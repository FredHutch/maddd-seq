#!/usr/bin/env python3

from collections import defaultdict
import gzip
import os
import pysam
import re

# Parse the params provided by Nextflow

# BARCODE LENGTH
barcode_length = int("${params.barcode_length}")
print(f"Removing {barcode_length} bases from the 5' of each read")

# MAXIMUM HOMOPOLYMER LENGTH
max_homopolymer = int("${params.barcode_max_homopolymer}")
print(f"Removing any barcodes with homopolymers longer than {max_homopolymer}bps")

# MINIMUM READ LENGTH
min_read_length = int("${params.min_align_score}")
print(f"Removing any reads which are shorter than {min_read_length} after barcode removal")


# The files are provided by Nextflow as R1 and R2
input_1 = "${R1}"
assert os.path.exists(input_1)
print(f"R1: {input_1}")
input_2 = "${R2}"
assert os.path.exists(input_2)
print(f"R2: {input_2}")

# Keep a counter of the input and output for both files
input_counter = 0
output_counter = 0

output_1 = "${R1.name.replaceAll(/.fastq.gz/, '')}.clipped.fastq.gz"
print(f"Output R1: {output_1}")
output_2 = "${R2.name.replaceAll(/.fastq.gz/, '')}.clipped.fastq.gz"
print(f"Output R2: {output_2}")

def has_homopolymer(s, n, chars=['A', 'T', 'C', 'G', 'a', 't', 'c', 'g']):
    """Check if a string s contains a homopolymer of length n"""

    # Iterate over the nucleotides
    for c in chars:

        # Check to see if there is a repeated character
        if re.search(r'%s{%s}' % (c, n), s) is not None:

            # If so, return True
            return True

    # If none of the characters matched, return False
    return False

def fails_quality_filter(bc):
    """Return True if the barcode sequence fails quality checks"""

    # If the barcode contains a homopolymer of length max_homopolymer + 1
    if has_homopolymer(bc, max_homopolymer + 1):

        # It fails
        return False


def reformat_read(read_x, bc_concat):
    """Reformat a single read"""

    # Add the concatenated barcode to the comment line
    read_x.comment = f"BC:Z:{bc_concat}"

    # Remove the barcode sequence
    read_x.sequence = read_x.sequence[barcode_length:]

    # Remove the barcode quality
    read_x.quality = read_x.quality[barcode_length:]

    return read_x


# Open all of the input and output file handles
with pysam.FastxFile(input_1) as i_1, pysam.FastxFile(input_2) as i_2, gzip.open(output_1, 'wt') as o_1, gzip.open(output_2, 'wt') as o_2:

    # Iterate over each of the reads in R1 and R2
    for read1, read2 in zip(i_1, i_2):

        # Add to the counter
        input_counter += 1

        # If we have processed 1M sequences
        if input_counter % 1000000 == 0:
            print(f"Processed {input_counter:,} read pairs")

        # Make sure that the reads are paired
        assert read1.name == read2.name, "Reads must be paired"

        # Get the first `barcode_length` bases from each read
        bc1 = read1.sequence[:barcode_length]
        bc2 = read2.sequence[:barcode_length]

        # If either barcode fails the quality filter
        if fails_quality_filter(bc1) or fails_quality_filter(bc2):
            
            # Skip it
            continue

        # Concatenate the barcode sequences
        bc_concat = bc1 + bc2

        # Reformat both reads, removing the barcode sequence and adding
        # the concatenated sequence to both
        read1 = reformat_read(read1, bc_concat)
        read2 = reformat_read(read2, bc_concat)

        # If either read is too short
        if len(read1.sequence) < min_read_length or len(read2.sequence) < min_read_length:

            # Skip it
            continue

        # Write out both reads
        # double slashes are used below to account for Nextflow interpolation
        o_1.write(str(read1) + "\\n")
        o_2.write(str(read2) + "\\n")

        output_counter += 1

print(f"Processed {input_counter:,} read pairs")
print(f"Wrote out {output_counter:,} read pairs passing all filters")
