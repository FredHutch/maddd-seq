#!/usr/bin/env python3
"""
Filter out any re-aligned sequences whose anchor (outermost) coordinates
have changed after the realignment of SSCs
"""

from collections import defaultdict
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
stats_output_fp = "${specimen}.unfiltered.SSC.details.csv.gz"
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

    # Populate a set of read names which have been seen
    seen = defaultdict(set)

    # Populate a set of read names to keep
    to_keep = defaultdict(set)

    # Populate a set of read names which are found multiple times
    multiples = set()

    # Keep a set of all of the index positions which have inconsistent positions
    inconsistent = set()

    # Initialize the counter
    i = 0
    
    # Open the BAM file for reading
    with pysam.AlignmentFile(bam_fp, "rb") as bam:

        # Iterate over each read
        for read in bam:

            # Add the read to the appropriate set based on its expected and observed position
            tally_read(read, to_keep, multiples, inconsistent, seen)

            # Increment the counter
            i += 1

    assert i > 0, "No reads found in input"

    print(f"Found {len(inconsistent):,} read pairs aligning to an unexpected position")
    print(f"Found {len(multiples):,} read pairs aligning multiple times")

    # Keep the reads which are aligned in the appropriate position
    # for both the forward and reverse reads
    to_keep = to_keep['fwd'].intersection(to_keep['rev'])

    # Combine the inconsistent and multiply-aligning reads and omit them all
    to_omit = inconsistent.union(multiples)
    print(f"Omitting {len(to_omit):,} read pairs overall")

    # Remove the inconsistent and multiple reads from the to_keep sets
    to_keep = to_keep - to_omit
    print(f"Found {len(to_keep):,} read pairs to keep")

    return to_keep, to_omit


def tally_read(read, to_keep, multiples, inconsistent, seen):
    """Add the read to the appropriate set based on its expected and observed position."""

    # Parse the barcode, chromosome, and position from the read header
    # Note that positions in the family ID are 1-indexed
    _, cname, fwd_pos, rev_pos = read.query_name.split("-")

    # Convert the family position information to 0-index, to conform
    # with how PySam parses the BAM file
    fwd_pos = int(fwd_pos) - 1
    rev_pos = int(rev_pos) - 1

    # Set up the key used to reference the dictionaries used to keep track of reads
    key = "rev" if read.is_reverse else "fwd"

    # If this read has been seen before
    if read.query_name in seen[key]:

        # Mark it as a multiple-aligner
        multiples.add(read.query_name)
        return

    # Mark this read as seen
    seen[key].add(read.query_name)

    # If the chromosome does not match the expectation
    if read.reference_name != cname:

        # Mark the read as inconsistent
        inconsistent.add(read.query_name)
        return

    # Calculate the offset from the expected position

    # If the read is in the reverse direction
    if read.is_reverse:

        # Compare the end of the alignment to the end of the family
        offset = abs(read.reference_end - rev_pos)

    # If the read is in the forward direction
    else:

        # If the start of the alignment matches the expectation
        offset = abs(read.reference_start - fwd_pos)

    # If the offset is below the threshold
    if offset <= max_realign_offset:

        # Add it to the set to keep
        to_keep[key].add(read.query_name)
    
    # Otherwise, if the offset is above the threshold
    else:

        # Add it to the set of inconsistent reads
        inconsistent.add(read.query_name)


def filter_bam(keep_reads, input_bam, output_bam):
    """Filter the input BAM to a set of reads."""

    print(f"Preparing to write {len(keep_reads):,} read pairs from {input_bam} to {output_bam}")

    # Start a counter
    i = 0

    # Open the input BAM file for reading
    with pysam.AlignmentFile(input_bam, "rb") as bam_i:
    
        # Open the output BAM file for writing
        with pysam.AlignmentFile(output_bam, "wb", template=bam_i) as bam_o:

            fwd_counter = 0
            rev_counter = 0

            # Iterate over the input
            for read in bam_i:

                # If the read is part of the set
                # and this position has not been masked due to inconsistent position
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

                # Increment the counter
                i += 1

    print(f"Wrote {fwd_counter:,} forward reads from {input_bam} to {output_bam}")
    assert len(keep_reads) == fwd_counter
    print(f"Wrote {rev_counter:,} reverse reads from {input_bam} to {output_bam}")
    assert len(keep_reads) == rev_counter


# First, read in the CSV data and join together
# This function will also return the CSV with all combined stats
ssc_stats = combine_ssc_stats()

# Get the set of families whose coordinates have stayed the same
# for each strand independently
keep_families_pos, omit_families_pos = find_consistent_families(input_pos_bam)
keep_families_neg, omit_families_neg = find_consistent_families(input_neg_bam)

# Combine the set of reads which were omitted on either strand
to_omit = omit_families_pos.union(omit_families_neg)

# Filter the final outputs to just those families which are consistent on one
# strand and not-inconsistent on the other (to allow us to keep families which
# only align to one of the two strands)
keep_families_pos = keep_families_pos - to_omit
keep_families_neg = keep_families_neg - to_omit

filter_bam(keep_families_pos, input_pos_bam, "POS.SSC.bam")
filter_bam(keep_families_neg, input_neg_bam, "NEG.SSC.bam")

# If the reads from one strand did not pass the realignment filter, then
# we need to update the ssc_stats table accordingly

ssc_stats.loc[
    # Filter to those families which did not pass the positive filter
    ~ssc_stats.family.isin(
        keep_families_pos
    ),
    # Update the fields for the positive strand
    [
        "R1-fwd-n", "R1-fwd-len", "R2-rev-n", "R2-rev-len"
    ]
# Set those values to 0
] = 0

ssc_stats.loc[
    # Filter to those families which did not pass the negative filter
    ~ssc_stats.family.isin(
        keep_families_neg
    ),
    # Update the fields for the negative strand
    [
        "R2-fwd-n", "R2-fwd-len", "R1-rev-n", "R1-rev-len"
    ]
# Set those values to 0
] = 0

# Filter the SSC statistics and write out in CSV format
ssc_stats = ssc_stats.query(
    "`R1-fwd-n` > 0 or `R2-fwd-n` > 0"
)
print(f"After filtering, the CSV contains {ssc_stats.shape[0]:,} families")

# Make sure that the table contains the expected number of rows
assert len(keep_families_pos | keep_families_neg) == ssc_stats.shape[0]

# Split up the family label into its components: barcode, refname, start_pos, end_pos
ssc_stats = pd.concat(
    [
        ssc_stats,
        ssc_stats.family.apply(
            lambda family_id: pd.Series(dict(zip(['barcode', 'refname', 'start_pos', 'end_pos'], family_id.split("-"))))
        )
    ],
    axis=1
)

# Change the column names to more informatively refer to the positive
# and negative strands
ssc_stats = ssc_stats.rename(
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
).sort_values(
    by=["refname", "start_pos", "end_pos"]
)

print(f"Writing out to {stats_output_fp}")
ssc_stats.to_csv(stats_output_fp, index=None)
print("DONE")
