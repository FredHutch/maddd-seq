#!/usr/bin/env python3
"""
Perform error correction on a table of barcodes

A barcode will be assigned to a 'corrected' sequence if:

- The 'corrected' sequence has a greater number of read pairs assigned
- The 'corrected' sequence is no more than ${params.max_barcode_mismatch} bases different
"""

import pandas as pd
import re

# MAXIMUM HOMOPOLYMER LENGTH
max_homopolymer = int("${params.barcode_max_homopolymer}")
print(f"Removing any barcodes with homopolymers longer than {max_homopolymer}bps")

# MAXIMUM NUMBER OF MISMATCHES BETWEEN BARCODES
max_barcode_mismatch = int("${params.max_barcode_mismatch}")
print(f"Maximum barcode mismatch: {max_barcode_mismatch}")

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


def has_homopolymer(s, n, chars=['A', 'T', 'C', 'G']):
    """Check if a string s contains a homopolymer of length n"""

    # Iterate over the nucleotides
    for c in chars:

        # Check to see if there is a repeated character
        if re.search(r'%s{%s}' % (c, n), s) is not None:

            # If so, return True
            return True

    # If none of the characters matched, return False
    return False


def passes_quality_filter(bc):
    """Return True if the barcode sequence passes quality checks"""

    # If the barcode contains any N's
    if "N" in bc or "n" in bc:

        # It fails
        return True

    # If the barcode contains a homopolymer of length max_homopolymer + 1
    if has_homopolymer(bc, max_homopolymer + 1):

        # It fails
        return True

    # If passes
    return False


# Read in the table of counts
print("Reading in barcode_counts.csv.gz")
df = pd.read_csv("barcode_counts.csv.gz")

print(f"Read in a table of {df.shape[0]:,} barcodes over {df['count'].sum():,} read pairs")

# Remove barcodes which fail the quality check
df = df.loc[df.barcode.apply(passes_quality_filter)]
print(f"After filtering on homopolymers and N's {df.shape[0]:,} barcodes remain over {df['count'].sum():,} read pairs")

# Sort by the number of read pairs
df.sort_values(by='count', inplace=True)

# Also make a copy of the table which is in desending order of counts
df_desc = df.sort_values(by='count', ascending=False)

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
