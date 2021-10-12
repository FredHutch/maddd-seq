#!/usr/bin/env python3
"""
Filter out any re-aligned sequences whose anchor (outermost) coordinates
have changed after the realignment of SSCs
"""

import os
import pandas as pd

# Input filepaths for realigned SSC sequences
input_pos_bam = "realigned.POS.SSC.bam"
assert os.path.exists(input_pos_bam)
input_neg_bam = "realigned.NEG.SSC.bam"
assert os.path.exists(input_neg_bam)


def combine_ssc_stats():
    """Function to read in and save SSC stats from all of the shards."""

    # Folder with input CSV data from all shards
    stats_folder = "SSC_STATS"
    assert os.path.exists(stats_folder)

    # List of files to read
    csv_list = [
        os.path.join(stats_folder, fp)
        for fp in os.listdir(stats_folder)
        if fp.endswith(".csv.gz")
    ]
    print(f"Reading in SSC stats from {len(csv_list):,} files")

    # Read in and combine all tables
    ssc_stats = pd.concat([pd.read_csv(fp) for fp in csv_list])

    # Output filepath
    stats_output_fp = "SSC.details.csv.gz"
    print(f"Writing out to {stats_output_fp}")

    # Write out in CSV format
    ssc_stats.to_csv(stats_output_fp, index=None)

    # Also return the DataFrame to the larger scope
    return ssc_stats

# First, read in the CSV data and join together
# This function will also return the CSV with all combined stats
ssc_stats = combine_ssc_stats()

