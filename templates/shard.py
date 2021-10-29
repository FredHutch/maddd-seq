#!/usr/bin/env python3

import gzip
import os
import pysam
import pandas as pd

read1_fp = '${R1}'
print(f"Read 1 FASTQ: {read1_fp}")
assert os.path.exists(read1_fp)

read2_fp = '${R2}'
print(f"Read 2 FASTQ: {read2_fp}")
assert os.path.exists(read2_fp)

# Read in the table of corrected barcodes
corrected_barcode_fp = "barcode_corrections.csv.gz"

corrected_barcodes = pd.read_csv(
    corrected_barcode_fp
).set_index(
    'barcode'
)[
    'corrected'
].to_dict()

print(f"Read in {len(corrected_barcodes):,} barcode corrections")

# Interpolate the number of shards from the nextflow param
# this will look odd when linting this particular file, but
# the syntax below will be replaced by nextflow prior to execution by Python
n_shards = ${params.n_shards}
print(f"Splitting into {n_shards:,} shards")

# Assign each corrected barcode to a shard index
shard_ix = dict()

# Build a set with the unique barcodes observed so far
unique_barcodes = set()

# Open both of the input files
with pysam.FastxFile(read1_fp) as i_1, pysam.FastxFile(read2_fp) as i_2:

    # Open a file handle for each of the shards
    print(f"Opening {n_shards:,} output file handles")
    output_handles = dict()

    # Iterate over 0..n_shards
    for i in range(n_shards):

        # Open the output file handles, one for each read
        output_handles[i] = [
            gzip.open(f"shard.{i}.R{R}.fastq.gz", "wt")
            for R in [1, 2]
        ]

    counter = 0

    # Iterate over the reads
    for read1, read2 in zip(i_1, i_2):

        # The ID for each read should match
        assert read1.name == read2.name, (read1.name, read2.name)

        # The barcode for each read should match
        assert read1.comment == read2.comment, (read1.comment, read2.comment)

        # Get the corrected barcode
        read_corrected_barcode = corrected_barcodes.get(read1.comment)

        # If this barcode could not be assigned to a corrected form
        if read_corrected_barcode is None:

            # Skip it
            continue

        # If we haven't seen this barcode before
        if read_corrected_barcode not in unique_barcodes:

            # Then we can assign a new index for the barcode
            # in the range(n_shards)
            shard_ix[read_corrected_barcode] = len(unique_barcodes) % n_shards

            # And add it to the set
            unique_barcodes.add(read_corrected_barcode)

        # For each of the reads
        for read_i, read in enumerate([read1, read2]):

            # Set the corrected barcode
            read.comment = read_corrected_barcode

            # Write out to the output handle
            output_handles[
                shard_ix[read_corrected_barcode]
            ][
                read_i
            ].write(
                str(read) + "\\n"
            )

            # Increment the counter
            counter += 1

print(f"Wrote out {counter:,} read pairs")
assert counter > 0

# Close all of the file handles
print(f"Closing {n_shards:,} output file handles")
for i, handle_list in output_handles.items():
    # Close both the R1 and the R2 output handles
    for handle in handle_list:
        handle.close()

print("DONE")