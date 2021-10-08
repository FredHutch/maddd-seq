#!/usr/bin/env python3

import gzip
import os
import pysam

bam_fp = '${bam}'
print(f"Input BAM: {bam_fp}")
assert os.path.exists(bam_fp)

barcodes_csv_gz = '${barcodes_csv_gz}'
print(f"Input CSV: {barcodes_csv_gz}")
assert os.path.exists(barcodes_csv_gz)

# Interpolate the number of shards from the nextflow param
# this will look odd when linting this particular file, but
# the syntax below will be replaced by nextflow prior to execution by Python
n_shards = ${params.n_shards}
print(f"Splitting into {n_shards:,} shards")

def yield_csv_gz(fp):
    """Yield each row of a CSV as a dict."""

    # Open the file
    with gzip.open(fp, 'rt') as handle:
    
        # Get the header
        header = handle.readline().rstrip("\\n").split(",")

        # Iterate over the subsequent lines
        for i, line in enumerate(handle):

            # Yield the row as a dict
            yield i, dict(zip(header, line.rstrip("\\n").split(",")))

# Assign each corrected barcode to a shard index
shard_ix = dict()

# Iterate over each row in the CSV
for i, r in yield_csv_gz(barcodes_csv_gz):

    # The corrected barcode sequence starts with "BC:Z"
    assert r["corrected"].startswith("BC:Z:")
    corrected = r["corrected"][len("BC:Z:"):]

    # If we haven't assigned this corrected barcode yet
    if shard_ix.get(corrected) is None:

        # Assign the shard IX based on the row number
        shard_ix[corrected] = i % n_shards

print(f"Assigned shards for {len(shard_ix)} barcodes")

# Open the input file
with pysam.AlignmentFile(bam_fp, "r") as bam:

    # Open a file handle for each of the shards
    print(f"Opening {n_shards:,} output file handles")
    output_handles = dict()

    # Iterate over 0..n_shards
    for i in range(n_shards):

        # Open the output file handle, using the header from the input
        output_handles[i] = pysam.AlignmentFile(
            f"shard.{i}.bam", 
            "wb",
            template=bam
        )

    counter = 0
    filtered = 0

    # Iterate over the reads
    for read in bam:

        # Get the shard for this barcode
        ix = shard_ix.get(read.get_tag('BC'))

        # If there is one
        if ix is not None:

            # Write out to the file for the shard
            output_handles[ix].write(read)
            counter += 1

        else:
            filtered += 1

print(f"Wrote out {counter:,} reads")
assert counter > 0
print(f"Filtered out {filtered:,} reads with non-matching barcodes")

# Close all of the file handles
print(f"Closing {n_shards:,} output file handles")
for i, handle in output_handles.items():
    handle.close()

print("DONE")