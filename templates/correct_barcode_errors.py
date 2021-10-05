#!/usr/bin/env python3
"""
Perform error correction on a table of barcodes

A barcode will be assigned to a 'corrected' sequence if:

- The 'corrected' sequence has a greater number of read pairs assigned
- The 'corrected' sequence is no more than ${params.max_barcode_mismatch} bases different
"""

import pandas as pd

max_barcode_mismatch = int("${params.max_barcode_mismatch}")
print(f"Maximum barcode mismatch: {max_barcode_mismatch}")

# Read in the table of counts
print("Reading in barcode_counts.csv.gz")
df = pd.read_csv("barcode_counts.csv.gz")

print(f"Read in a table of {df.shape[0]:,} barcodes over {df['count'].sum():,} read pairs")

# Sort by the number of read pairs
df.sort_values(by='count', inplace=True)

# Also make a copy of the table which is in desending order of counts
df_desc = df.sort_values(by='count', ascending=False)

def hamming_distance(s1, s2):
    """Return the hamming distance, the number of characters which do not match between two strings."""

    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


def find_corrected_barcode(bc_string, n):
    """Find a corrected barcode based on hamming distance and counts."""

    # Iterate over all of the barcodes which have greater counts
    for _, r in df_desc.query(f"count > {n}").iterrows():

        # If the sequences are sufficiently similar
        if hamming_distance(bc_string, r["barcode"]) <= max_barcode_mismatch:

            # Then return the corrected sequence
            return r["barcode"]

    # If no matches were found

    # Return the original sequence
    return r["barcode"]


# Make a dictionary linking each barcode to its error-corrected form
corrected = dict()

# Iterate over each barcode, moving from the lowest-abundance to the greatest
for _, r in df.iterrows():

    # Assign the corrected barcode, if any
    corrected[r['barcode']] = find_corrected_barcode(r['barcode'], r['count'])

# Add the corrected barcode to the table
df = df.assign(
    corrected=df['barcode'].apply(corrected.get)
).sort_values(
    by="count",
    ascending=False
)

# Write to the file
df.to_csv("barcode_corrections.csv.gz", index=None)
