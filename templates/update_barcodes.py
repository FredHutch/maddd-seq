#!/usr/bin/env python3
"""Apply the corrected barcode sequences to the WGS data in FASTQ format"""

import gzip
import os
import pysam

# The files are provided by Nextflow as R1 and R2
input_1 = "${R1}"
assert os.path.exists(input_1)
print(f"R1: {input_1}")
input_2 = "${R2}"
assert os.path.exists(input_2)
print(f"R2: {input_2}")

output_1 = "${R1.name.replaceAll(/.fastq.gz/, '')}.corrected.fastq.gz"
print(f"Output R1: {output_1}")
output_2 = "${R2.name.replaceAll(/.fastq.gz/, '')}.corrected.fastq.gz"
print(f"Output R2: {output_2}")


# Read in the table of corrected barcodes
corrected_barcode_fp = "barcode_corrections.csv.gz"

# Get the header
with gzip.open(corrected_barcode_fp, "rt") as handle:
    table_header = handle.readlines(1).rstrip("\\n").split(",")

# Make sure that we have the headers that are expected, and get the index
assert 'barcode' in table_header
assert 'corrected' in table_header

# Make a dict with the corrected sequences
corrected_barcodes = dict()

# Iterate over the table
with gzip.open(corrected_barcode_fp, "rt") as handle:
    for i, l in enumerate(handle):

        # Skip the first line
        if i == 0:
            continue

        # Parse the line
        fields = dict(
            zip(
                table_header,
                l.rstrip("\\n").split(",")
            )
        )

        # Add the data from this row to the dict
        corrected_barcodes[
            fields['barcode']
        ] = fields['corrected']

def reformat_read(read_x):
    """Reformat a single read"""

    # Get the original barcode sequence
    assert read_x.comment.startswith("BC:Z:")
    original_bc = read_x.comment[len("BC:Z:"):]

    # Get the corrected barcode
    assert original_bc in corrected_barcodes
    new_bc = corrected_barcodes[original_bc]

    # Fill in the corrected barcode
    read_x.comment = f"BC:Z:{new_bc}"

    return read_x


def correct_fastq(fp_in, fp_out):
    """Correct all of the barcodes in a FASTQ file."""

    counter = 0

    # Open all of the input and output file handles
    with pysam.FastxFile(fp_in) as handle_i, gzip.open(fp_out, 'wt') as handle_o:

        # Iterate over each of the reads
        for read in handle_i:

            # Correct the read and write out
            handle_o.write(
                str(
                    reformat_read(read)
                ) + "\\n"
            )
            counter += 1

    print(f"Corrected {counter:,} reads from {fp_in} --> {fp_out}")

# Process both R1 and R2 in the same way
correct_fastq(input_1, output_1)
correct_fastq(input_2, output_2)
