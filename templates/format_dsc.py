#!/usr/bin/env python3
"""
Parse the SSC data in order to generate a BAM with the double-strand consensus (DSC)
"""

import gzip
import json
import os
import pandas as pd
import pysam

# Input filepaths for filtered SSC sequences
input_pos_bam = "POS.SSC.bam"
assert os.path.exists(input_pos_bam)
input_neg_bam = "NEG.SSC.bam"
assert os.path.exists(input_neg_bam)

# Output path for the DSC BAM
dsc_output_unsorted = "unsorted.DSC.bam"
dsc_output = "DSC.bam"
print(f"Output path: {dsc_output}")


# Define the rules for combining two nucleotides
# M = A or C
# R = A or G
# W = A or T
# S = C or G
# Y = C or T
# K = G or T
iupac = dict(
    A=dict(A='A', T='W', C='M', G='R', N='N'),
    T=dict(A='W', T='T', C='Y', G='K', N='N'),
    C=dict(A='M', T='Y', C='C', G='S', N='N'),
    G=dict(A='R', T='K', C='S', G='G', N='N'),
    N=dict(A='N', T='N', C='N', G='N', N='N')
)

def combine_bases(base_a, base_b):
    """Return the appropriate nucleotide based on the two bases read from each strand."""
    return iupac[base_a.upper()][base_b.upper()]

def compute_consensus(pos_ssc, neg_ssc):
    """Compute the consensus nucleotide sequence from two sequences."""

    # Make sure that the read IDs match
    assert pos_ssc.query_name == neg_ssc.query_name, (pos_ssc.query_name, neg_ssc.query_name)

    # Replace the sequence of one read with the consensus
    pos_ssc.query_sequence = ''.join([
        combine_bases(base_a, base_b)
        for base_a, base_b in zip(pos_ssc.query_sequence, neg_ssc.query_sequence)
    ])

    # Return the consensus read
    return pos_ssc

# Open both BAM files for the input
with pysam.AlignmentFile(input_pos_bam, "rb") as pos_bam, pysam.AlignmentFile(input_neg_bam, "rb") as neg_bam:

    # Open the output, using the header from the input
    # Write to the unsorted file path, which will be sorted later
    with pysam.AlignmentFile(dsc_output_unsorted, "w", template=pos_bam) as dsc_bam:

        # Iterate over each set of inputs
        for pos_ssc, neg_ssc in zip(pos_bam, neg_bam):

            # Write out the combined output
            dsc_bam.write(
                compute_consensus(pos_ssc, neg_ssc)
            )

# Now sort the output
pysam.sort("-o", dsc_output, dsc_output_unsorted)

# And finally index it
pysam.index(dsc_output)

print("DONE")