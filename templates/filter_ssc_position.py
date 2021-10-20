#!/usr/bin/env python3
"""
Filter out any re-aligned sequences whose anchor (outermost) coordinates
have changed after the realignment of SSCs
"""

import os
import pandas as pd
import pysam
from pandas.errors import EmptyDataError

# Input filepaths for realigned SSC sequences
input_pos_bam = "realigned.POS.SSC.bam"
assert os.path.exists(input_pos_bam)
input_neg_bam = "realigned.NEG.SSC.bam"
assert os.path.exists(input_neg_bam)

# Output filepath for SSC details
stats_output_fp = "SSC.details.csv.gz"
print(f"Writing out to {stats_output_fp}")

# Maximum offset following realignment
max_realign_offset = int("${params.max_realign_offset}")
print(f"Maximum offset after realignment: {max_realign_offset}bp")


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
    ssc_stats = try_read_csv_list(csv_list)

    # Return the DataFrame
    return ssc_stats


def try_read_csv_list(csv_list):
    """Try to read and concatenate a list of CSV filepaths."""

    # Populate a list of DataFrames
    df_list = list()

    # Iterate over the filepaths
    for fp in csv_list:

        # Try to read the file
        try:
            df = pd.read_csv(fp)

        # If there is no data
        except EmptyDataError:

            # Skip the file
            continue

        # If there was no exception
        # Add the data to the list
        df_list.append(df)

    # Concatenate all of the data
    return pd.concat(df_list)


def find_consistent_families(bam_fp):
    """Return the set of families whose positions have stayed consistent after re-mapping."""

    # Populate a set of read names to keep
    to_keep = set()

    # Initialize the counter
    i = 0
    
    # Open the BAM file for reading
    with pysam.AlignmentFile(bam_fp, "rb") as bam:

        # Iterate over each read
        for read in bam:

            # Increment the counter
            i += 1

            # Parse the barcode, chromosome, and position from the read header
            _, cname, fwd_pos, rev_pos = read.query_name.split("-")

            # If the chromosome matches the expectation
            if read.reference_name == cname:

                # If the read is in the reverse direction
                if read.is_reverse:

                    # If the end of the alignment matches the expectation
                    if abs(read.reference_end - int(rev_pos)) <= max_realign_offset:

                        # Add it to the set to keep
                        to_keep.add(read.query_name)

                    else:
                        print(read.query_name)
                        print('rev' if read.is_reverse else 'fwd')
                        print(read.reference_start)
                        print(read.reference_end)

                # If the read is in the forward direction
                else:

                    # If the start of the alignment matches the expectation
                    if abs(read.reference_start - int(rev_pos)) <= max_realign_offset:

                        # Add it to the set to keep
                        to_keep.add(read.query_name)

                    else:
                        print(read.query_name)
                        print('rev' if read.is_reverse else 'fwd')
                        print(read.reference_start)
                        print(read.reference_end)

    assert i > 0, "No reads found in input"

    print(f"Keeping {len(to_keep):,} / {i+1:,} reads from {bam_fp}")
    return to_keep


def filter_bam(keep_reads, input_bam, output_bam):
    """Filter the input BAM to a set of reads."""

    print(f"Preparing to write {len(keep_reads):,} read pairs from {input_bam} to {output_bam}")

    # Open the input BAM file for reading
    with pysam.AlignmentFile(input_bam, "rb") as bam_i:
    
        # Open the output BAM file for writing
        with pysam.AlignmentFile(output_bam, "wb", template=bam_i) as bam_o:

            fwd_counter = 0
            rev_counter = 0

            # Iterate over the input
            for read in bam_i:

                # If the read is part of the set
                if read.query_name in keep_reads:

                    # Write it out
                    bam_o.write(read)

                    # If the read is reverse
                    if read.is_reverse:

                        # Increment the reverse counter
                        rev_counter += 1

                    # If the read is forward
                    else:

                        # Increment the forward counter
                        fwd_counter += 1

    print(f"Wrote {fwd_counter:,} forward reads from {input_bam} to {output_bam}")
    assert len(keep_reads) == fwd_counter
    print(f"Wrote {rev_counter:,} reverse reads from {input_bam} to {output_bam}")
    assert len(keep_reads) == rev_counter


# First, read in the CSV data and join together
# This function will also return the CSV with all combined stats
ssc_stats = combine_ssc_stats()

# Get the set of families whose coordinates have stayed the same
# for each strand independently
keep_families_pos = find_consistent_families(input_pos_bam)
keep_families_neg = find_consistent_families(input_neg_bam)

# Find the intersection, families which are consistent on both strands
keep_families = keep_families_pos & keep_families_neg
print(f"Keeping {len(keep_families):,} reads from both strands")

# Filter the final outputs to just that set of families
filter_bam(keep_families, input_pos_bam, "POS.SSC.bam")
filter_bam(keep_families, input_neg_bam, "NEG.SSC.bam")

# Filter the SSC statistics and write out in CSV format
ssc_stats = ssc_stats.loc[
    ssc_stats.family.isin(keep_families)
]
print(f"Writing out to {stats_output_fp}")
ssc_stats.to_csv(stats_output_fp, index=None)
print("DONE")
