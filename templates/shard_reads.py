#!/usr/bin/env python3

import gzip
import os
import pysam
import pandas as pd

specimen = '${specimen}'
print(f"Specimen: {specimen}")

shard_ix = '${shard_ix}'
print(f"Shard: {shard_ix}")

read1_fp = '${R1}'
print(f"Read 1 FASTQ: {read1_fp}")
assert os.path.exists(read1_fp)

read2_fp = '${R2}'
print(f"Read 2 FASTQ: {read2_fp}")
assert os.path.exists(read2_fp)

barcode_fp = '${shard_barcodes_csv}'
print(f"Barcodes: {barcode_fp}")
assert os.path.exists(barcode_fp)

# Read in the table of corrected barcodes
corrected_barcodes = pd.read_csv(
    barcode_fp
).set_index(
    'barcode'
)[
    'corrected'
].to_dict()

print(f"Read in {len(corrected_barcodes):,} barcode corrections")

# Open both of the input files
with pysam.FastxFile(read1_fp) as i_1, pysam.FastxFile(read2_fp) as i_2:

    # Open the output file handles, one for each read
    output_handles = [
        gzip.open(f"{shard_ix}.R{R}.fastq.gz", "wt")
        for R in [1, 2]
    ]

    counter = 0

    # Iterate over the reads
    for read1, read2 in zip(i_1, i_2):

        # Get the corrected barcode
        read_corrected_barcode = corrected_barcodes.get(read1.comment)

        # If this barcode could not be assigned to a corrected form
        if read_corrected_barcode is None:

            # Skip it
            continue

        # The ID for each read should match
        assert read1.name == read2.name, (read1.name, read2.name)

        # The barcode for each read should match
        assert read1.comment == read2.comment, (read1.comment, read2.comment)

        # For each of the reads
        for read_i, read in enumerate([read1, read2]):

            # Set the corrected barcode
            read.comment = read_corrected_barcode

            # Write out to the output handle
            output_handles[
                read_i
            ].write(
                str(read) + "\\n"
            )

            # Increment the counter
            counter += 1

print(f"Wrote out {counter:,} read pairs")
assert counter > 0

# Close both the R1 and the R2 output handles
for handle in output_handles:
    handle.close()

print("DONE")