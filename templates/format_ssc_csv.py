#!/usr/bin/env python3

import gzip
import json
import pandas as pd

# Single JSON with all of the information for each DSC
with gzip.open('total.json.gz', 'rt') as handle:
    dat = json.load(handle)

# Get information for each DSC, including:
# - the merged length of forward and reverse reads
merged_len = pd.Series({k: v['nbases'] for k, v in dat.items()})
# - the number of adducts
n_adducts = pd.Series({k: len(v['adducts']) for k, v in dat.items()})
# - the number of mutations
n_mutations = pd.Series({k: len(v['mutations']) for k, v in dat.items()})

# Read the table with the number of reads per SSC
ssc_df = pd.read_csv("SSC.details.csv.gz")

# Get the number of reads on the positive and negative strands
# Get the length of the forward and reverse reads
ssc_df = ssc_df.rename(
    columns={
        "R1-fwd-n": "nreads_pos",
        "R1-fwd-len": "rlen_fwd",
        "R2-fwd-n": "nreads_neg",
        "R1-rev-len": "rlen_rev"
    }
).drop(
    columns=[
        f"{R}-{d}-{m}"
        for R in ['R1', 'R2']
        for d in ['fwd', 'rev']
        for m in ['n', 'len']
    ],
    errors='ignore'
).set_index(
    'family'
)

# Add the data on merged_len, n_adducts, and n_mutations
ssc_df = ssc_df.assign(
    merged_len=merged_len,
    n_adducts=n_adducts,
    n_mutations=n_mutations,
).reset_index()

# Split up the family name into barcode, chromosome, and positions
ssc_df = pd.concat(
    [
        pd.DataFrame([
            dict(zip(['barcode', 'chr', 'start', 'stop'], family_str.split("-")))
            for family_str in ssc_df['family'].values
        ]),
        ssc_df,
    ],
    axis=1
).drop(columns='family')

# Save to disk
ssc_df.to_csv("${specimen}.SSC.csv.gz", index=None)