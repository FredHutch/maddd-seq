#!/usr/bin/env python3
"""
Filter out any re-aligned sequences whose anchor (outermost) coordinates
have changed after the realignment of SSCs
"""

import os
import pandas as pd
import pysam

# Minimum number of reads per SSC
min_reads = int("${params.min_reads}")
print(f"Filtering to >={min_reads} reads per SSC on both strands")

# Input filepaths for unfiltered SSC sequences
input_pos_bam = "unfiltered.POS.SSC.bam"
assert os.path.exists(input_pos_bam)
input_neg_bam = "unfiltered.NEG.SSC.bam"
assert os.path.exists(input_neg_bam)

# Output filepaths for filtered SSC sequences
output_pos_bam = "POS.SSC.bam"
output_neg_bam = "NEG.SSC.bam"

# Input filepath for the SSC stats CSV
input_stats_csv = "unfiltered.SSC.details.csv.gz"
assert os.path.exists(input_stats_csv)

# Output filepath for filtered stats CSV
output_stats_csv = "SSC.details.csv.gz"

# Read in the stats with the number of reads for each part of each family
# In this case it's important to remember that each family has a R1/R2 and fwd/rev SSC
ssc_stats = pd.read_csv(input_stats_csv)
print(f"Read in stats for {ssc_stats.shape[0]:,} SSC families")


def filter_sscs(ssc_stats, min_reads):
    """
    Filter the SSCs based on the number of reads in both strands.
    The output will be two sets (fwd_reads, rev_reads), which
    contain the set of family IDs for those SSCs which should be
    kept in the FWD and REV directions, respectively.
    """

    fwd_reads = set(
        ssc_stats.loc[
            (ssc_stats["R1-fwd-n"] >= min_reads & ssc_stats["R2-fwd-n"] >= min_reads),
            'family'
        ].tolist()
    )
    rev_reads = set(
        ssc_stats.loc[
            (ssc_stats["R1-rev-n"] >= min_reads & ssc_stats["R2-rev-n"] >= min_reads),
            'family'
        ].tolist()
    )

    return fwd_reads, rev_reads


def filter_bam(input_bam, output_bam, fwd_reads, rev_reads):
    """Filter a BAM to just the reads in the defined sets (for forward and reverse reads independently)"""

    print(f"Preparing to write {len(fwd_reads):,} forward reads and {len(rev_reads):,} reverse reads")
    print(f"Writing from {input_bam} to {output_bam}")

    # Open the input BAM file for reading
    with pysam.AlignmentFile(input_bam, "rb") as bam_i:
    
        # Open the output BAM file for writing
        with pysam.AlignmentFile(output_bam, "wb", template=bam_i) as bam_o:

            fwd_counter = 0
            rev_counter = 0

            # Iterate over the input
            for read in bam_i:

                # If the read is oriented in the reverse direction
                if read.is_reverse:

                    # And it is in the set of reads to keep
                    if read.query_name in rev_reads:

                        # Write it out
                        bam_o.write(read)
                        rev_counter += 1

                # If it is oriented in the forward direction
                else:

                    # And it is in the set of reads to keep
                    if read.query_name in fwd_reads:

                        # Write it out
                        bam_o.write(read)
                        fwd_counter += 1

    print(f"Wrote out {fwd_counter:,} forward reads")
    assert len(fwd_reads) == fwd_counter
    print(f"Wrote out {rev_counter:,} reverse reads")
    assert len(rev_reads) == rev_counter


# Filter the SSCs based on the minimum number of reads on both strands
fwd_reads, rev_reads = filter_sscs(ssc_stats, min_reads)

# Filter the BAM files for both strands
filter_bam(input_pos_bam, output_pos_bam, fwd_reads, rev_reads)
filter_bam(input_neg_bam, output_neg_bam, fwd_reads, rev_reads)

# Filter the CSV
ssc_stats = ssc_stats.loc[
    ssc_stats['family'].isin(fwd_reads | rev_reads)
]
print(f"Writing out metrics for {ssc_stats.shape[0]:,} families to {output_stats_csv}")
assert ssc_stats.shape[0] == len(fwd_reads | rev_reads)
ssc_stats.to_csv(output_stats_csv, index=None)
