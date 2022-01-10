#!/usr/bin/env python3

import os
import pandas as pd

# Read in the table of corrected barcodes
corrected_barcode_fp = "barcode_corrections.csv.gz"
assert os.path.exists(corrected_barcode_fp)
barcode_df = pd.read_csv(corrected_barcode_fp)

# Sort the corrected barcodes by the number of reads assigned
corrected_barcodes = barcode_df.groupby(
    'corrected'
)[
    'count'
].sum(
).sort_values(
).index.values

print(f"Read in {len(corrected_barcodes):,} barcodes")

# Interpolate the number of shards from the nextflow param
# this will look odd when linting this particular file, but
# the syntax below will be replaced by nextflow prior to execution by Python
n_shards = ${params.n_shards}
print(f"Splitting into {n_shards:,} shards")

# Assign each corrected barcode to a shard
shard_dict = dict()
for barcode_i, barcode in enumerate(corrected_barcodes):
    shard_dict[barcode] = barcode_i % n_shards

# Add the shard to the barcode table
barcode_df = barcode_df.assign(
    shard = barcode_df['corrected'].apply(
        shard_dict.get
    )
)

# For each shard of barcodes
for shard_i, shard_barcodes in barcode_df.groupby('shard'):

    fp_out = f"shard.{shard_i}.barcode_corrections.csv.gz"

    print(f"Writing out {shard_barcodes.corrected.unique().shape[0]:,} barcodes to {fp_out}")

    # Write out the list of barcodes
    shard_barcodes.drop(
        columns=['shard']
    ).to_csv(
        fp_out,
        index=None
    )

print("DONE")
