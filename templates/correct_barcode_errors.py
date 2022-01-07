#!/usr/bin/env python3
"""
Perform error correction on a table of barcodes

A barcode will be assigned to a 'corrected' sequence if:

- The 'corrected' sequence has a greater number of read pairs assigned
- The 'corrected' sequence is no more than ${params.max_barcode_mismatch} bases different
"""

from functools import lru_cache
import logging
import os
import pandas as pd
from Bio.Seq import reverse_complement

# Set the level of the logger to INFO
logFormatter = logging.Formatter(
    '%(asctime)s %(levelname)-8s [correct_barcode_errors] %(message)s'
)
logger = logging.getLogger('correct_barcode_errors')
logger.setLevel(logging.INFO)

# Write to STDOUT
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
logger.addHandler(consoleHandler)

logger.info("Processing specimen ${specimen}")

# FILE CONTAINING BARCODES
barcodes_txt = "barcodes.txt"
logger.info(f"Reading barcodes from {barcodes_txt}")
assert os.path.exists(barcodes_txt)

# EXPECTED LENGTH OF BARCODES
barcode_len = ${params.barcode_length}
logger.info(f"All barcodes should be {barcode_len}bp")
assert isinstance(barcode_len, int)

# MAXIMUM NUMBER OF MISMATCHES BETWEEN BARCODES
max_barcode_mismatch = int("${params.max_barcode_mismatch}")
logger.info(f"Maximum barcode mismatch: {max_barcode_mismatch}")

fp_out = "barcode_corrections.csv.gz"
logger.info(f"Using output path: {fp_out}")


def read_barcodes(fp):
    """Read in the list of whitelist barcodes."""

    # Make a list of barcodes
    barcodes = list()

    # Open the file
    with open(fp, "r") as handle:

        # Iterate over each line
        for line in handle:

            # Strip off the newline character
            bc = line.rstrip("\\n")

            # If the line is empty
            if len(bc) == 0:

                # Skip it
                continue

            # Make sure that the barcode is the expected length
            assert len(bc) == barcode_len, len(bc)

            # Add it to the list
            barcodes.append(bc)

    return barcodes


def hamming_distance(s1, s2):
    """Return the hamming distance, the number of characters which do not match between two strings."""

    return min([
        sum(c1 != c2 for c1, c2 in zip(query, s2))
        for query in [s1, reverse_complement(s1)]
    ])


def correct_merged_barcode(bc_tag, tag_prefix="BC:Z:"):
    """Compute the corrected barcode for a merged barcode."""

    # The input barcode sequence starts with a prefix, which
    # should be removed prior to parsing
    assert bc_tag.startswith(tag_prefix)
    bc = bc_tag[len(tag_prefix):]

    # The length of the barcode should be 2x the individual barcode length
    assert len(bc) == barcode_len * 2, bc

    # Split up the barcode into each half, starting with the left half
    left_bc = correct_barcode(bc[:barcode_len])

    # If there was no match found
    if left_bc is None:

        # Then the merged barcode will be thrown out
        return None

    # Make sure that the corrected barcode is the correct length
    assert len(left_bc) == barcode_len, (bc[:barcode_len], left_bc)

    # Now correct the right half
    right_bc = correct_barcode(bc[barcode_len:])

    # If there was no match found
    if right_bc is None:

        # Then the merged barcode will be thrown out
        return None

    # Make sure that the corrected barcode is the correct length
    assert len(right_bc) == barcode_len, (bc[barcode_len:], right_bc)

    # At this point, there was a match for both
    # Return the merged barcode
    return tag_prefix + left_bc + right_bc


# Read in the table of counts
logger.info("Reading in barcode_counts.csv.gz")
df = pd.read_csv("barcode_counts.csv.gz")

logger.info(f"Read in a table of {df.shape[0]:,} barcodes over {df['count'].sum():,} read pairs")

# Read the list of barcodes
whitelist = read_barcodes(barcodes_txt)

logger.info(f"Read in a list of {len(whitelist):,} whitelist barcodes")


# The function below is being defined after reading in the whitelist
# so that the caching decorator can be used most effectively
@lru_cache
def correct_barcode(bc):

    # First check for an exact match
    for whitelist_bc in whitelist:
        if bc == whitelist_bc:
            return whitelist_bc

    # Next, calculate the hamming distance for each
    whitelist_distances = {
        whitelist_bc: hamming_distance(whitelist_bc, bc)
        for whitelist_bc in whitelist
    }

    # Get the lowest hamming distance
    lowest_distance = min(whitelist_distances.values())

    # If that distance meets the threshold
    if lowest_distance <= max_barcode_mismatch:

        # Get the whitelist barcodes which are at that minimum
        best_whitelist_bcs = [
            bc
            for bc, d in whitelist_distances.items()
            if d == lowest_distance
        ]

        # If there is only a single best match
        if len(best_whitelist_bcs) == 1:

            # Then return it
            return best_whitelist_bcs[0]

    # If there was no single best match which met the threshold
    return None


# Add a column to the DataFrame with the corrected barcode
# If there is no match to the whitelist, the populated value will be None
df = df.assign(
    corrected=df['barcode'].apply(
        lambda bc: correct_merged_barcode(bc)
    )
)

# Drop any barcodes which could not be corrected
df = df.dropna()
assert df.shape[0] > 0, "Could not find any barcode matches"
logger.info(f"Output: {df.shape[0]:,} barcodes with {df['count'].sum():,} total counts")

# Write to the file
logger.info(f"Writing out to {fp_out}")
df.to_csv(
    fp_out, 
    index=None
)

logger.info("DONE")
