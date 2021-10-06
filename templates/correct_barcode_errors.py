#!/usr/bin/env python3
"""
Perform error correction on a table of barcodes

A barcode will be assigned to a 'corrected' sequence if:

- The 'corrected' sequence has a greater number of read pairs assigned
- The 'corrected' sequence is no more than ${params.max_barcode_mismatch} bases different
"""

from collections import defaultdict
from bitarray import bitarray
import logging
import pandas as pd
import re

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

# MAXIMUM HOMOPOLYMER LENGTH
max_homopolymer = int("${params.barcode_max_homopolymer}")
logger.info(f"Removing any barcodes with homopolymers longer than {max_homopolymer}bps")

# MAXIMUM NUMBER OF MISMATCHES BETWEEN BARCODES
max_barcode_mismatch = int("${params.max_barcode_mismatch}")
logger.info(f"Maximum barcode mismatch: {max_barcode_mismatch}")

fp_out = "barcode_corrections.csv.gz"
logger.info(f"Using output path: {fp_out}")


def hamming_distance(s1, s2):
    """Return the hamming distance, the number of characters which do not match between two strings."""

    return sum(c1 != c2 for c1, c2 in zip(s1, s2))


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
        return False

    # If the barcode contains a homopolymer of length max_homopolymer + 1
    if has_homopolymer(bc, max_homopolymer + 1):

        # It fails
        return False

    # If passes
    return True


def make_barcode_index(barcodes):
    """Make a k-mer based index for a Series of barcodes."""

    # Get the length of the barcodes
    barcode_lengths = barcodes.apply(len)

    # All of the barcodes should be the same length
    barcode_lengths = barcode_lengths.unique()
    assert barcode_lengths.shape[0] == 1, "All barcodes must be the same length"
    barcode_length = barcode_lengths[0]
    logger.info(f"All barcodes are {barcode_length}bp")

    # Set the index size based on the number of mismatches
    index_k = int(barcode_length / max_barcode_mismatch)

    # The smallest value of k that we will allow is 4
    assert index_k >= 4, f"Cannot allow {max_barcode_mismatch} mismatches for barcodes of {barcode_length}bp"

    logger.info(f"Using an index size of {index_k}bp")

    # Transform the barcodes into a list of sets
    return [
        set(list(get_kmers(bc, index_k)))
        for bc in barcodes.values
    ]


def get_kmers(bc, k, step=4):
    "Yield the kmers of size `k` from string `bc`"
    
    for i in range(0, len(bc) - k, step):
        yield bc[i: (i + k)]


# Read in the table of counts
logger.info("Reading in barcode_counts.csv.gz")
df = pd.read_csv("barcode_counts.csv.gz")

logger.info(f"Read in a table of {df.shape[0]:,} barcodes over {df['count'].sum():,} read pairs")

# Remove barcodes which fail the quality check
df = df.loc[df.barcode.apply(passes_quality_filter)]
logger.info(f"After filtering on homopolymers and N's {df.shape[0]:,} barcodes remain over {df['count'].sum():,} read pairs")

# Sort by the number of read pairs
df = df.sort_values(by='count').reset_index(drop=True)
logger.info("Sorted barcodes by increasing counts")

# For rapid hashing, pull out the even and the odd positions
# and add them to a dedicated columns, 
# since we can only compare those barcodes which share
# a large amount of sequence (for compute limitations)
df = df.assign(
    evens=df.barcode.apply(lambda s: s[::2]),
    odds=df.barcode.apply(lambda s: s[1::2]),
)
logger.info("Annotated barcode prefixes and suffixes")

# Make a dictionary linking each barcode to its error-corrected form
corrected = dict()

# Join barcodes first by prefix, then by suffix
for cname in ['evens', 'odds']:
    logger.info(f"Joining barcodes by {cname}")

    # Get the counts for each prefix/suffix
    cname_vc = df[cname].value_counts()

    # Get the mask for barcodes with a prefix/suffix which is found >1 times
    cname_ix = df[cname].apply(cname_vc.get) > 1

    logger.info(f"Found {cname_ix.sum():,} barcodes with a shared {cname}")

    group_counter = 0
    bc_counter = 0

    # Iterate over each group
    for shared_seq, shared_df in df.loc[cname_ix].groupby(cname):

        group_counter += 1

        # For each barcode in this group
        for i, i_r in shared_df.iterrows():

            bc_counter += 1

            # If this barcode has already been corrected
            if i_r.barcode in corrected:

                # Skip it
                continue

            # If another barcode has been corrected to this one
            if i_r.barcode in corrected.values():

                # Skip it
                continue

            # If this barcode has not yet been corrected
            else:

                # Iterate over the other barcodes, starting from the greatest counts
                for j, j_r in shared_df[::-1].iterrows():

                    # If the other barcode has been corrected already
                    if j_r['barcode'] in corrected:

                        # Skip it
                        continue

                    # If the second barcode has greater counts
                    if j_r['count'] > i_r['count']:

                        # If the hamming distance is below the threshold
                        if hamming_distance(i_r['barcode'], j_r['barcode']) <= max_barcode_mismatch:

                            # The i'th barcode should be corrected to the j'th
                            corrected[i_r['barcode']] = j_r['barcode']

            if bc_counter > 0 and bc_counter % 10000 == 0:
                logger.info(f"Compared {bc_counter:,} barcodes across {group_counter:,} groups, assigned corrected barcodes for {len(corrected):,} barcodes")

logger.info(f"Assigned corrected barcodes for {len(corrected):,} barcodes")

# Add the corrected barcode to the table
df = df.assign(
    corrected=df['barcode'].apply(
        lambda s: corrected.get(s, s)
    )
).drop(
    columns=['evens', 'odds']
).sort_values(
    by="count",
    ascending=False
)

# Write to the file
logger.info(f"Writing out to {fp_out}")
df.to_csv(
    fp_out, 
    index=None
)

logger.info("DONE")
