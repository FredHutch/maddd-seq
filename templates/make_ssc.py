#!/usr/bin/env python3

from collections import defaultdict
from Bio import AlignIO
from Bio.Align.AlignInfo import SummaryInfo
from io import StringIO
import gzip
from math import ceil
import os
from Bio.Seq import reverse_complement
import pandas as pd
import pysam

# Parse the input filepaths

# CSV assigning each read to a family
families_csv_gz = "families.csv.gz"
print(f"Reading family assignments from {families_csv_gz}")
assert os.path.exists(families_csv_gz)

# Aligned reads in BAM format
bam = "${bam}"
print(f"bam: {bam}")
assert os.path.exists(bam)

# Minimum proportion of bases needed to calculate a consensus
min_base_prop = ${params.min_base_prop}
print(f"Minimum proportion of bases required for consensus: {min_base_prop}")

# Get the shard index, used to name the outputs
shard_ix = "${shard_ix}"

# Output file for the stats
stats_csv = f"{shard_ix}.stats.csv.gz"


# READ IN THE FAMILY ASSIGNMENT PER READ
def read_family_assignments(families_csv_gz):

    # Read in the DataFrame of families (id, family)
    families_df = pd.read_csv(families_csv_gz)
    print(f"Read in {families_df.shape[0]:,} lines from {families_csv_gz}")

    # Return the dictionary linking id -> family
    families = families_df.set_index('id')['family']
    # as well as the number of entries for each family
    return families, families.value_counts()


# PARSE THE INPUT READS - YIELD SSCs
def parse_input_bam(bam, families):

    # Keep the dicts, keyed by (1) the family ID and (2) R1/R2+fwd/rev
    family_reads = defaultdict(lambda: defaultdict(list))

    counter = 0

    # Open the BAM input file
    with pysam.AlignmentFile(bam, "rb") as bam:

        # Iterate over the reads
        for read in bam:

            # If this read did not have a family assigned,
            # due to a lack of unique pairing
            if families.get(read.query_name) is None:

                # Skip it
                continue

            # Otherwise, format the label by R1/R2 and fwd/rev
            side = "R1" if read.is_read1 else "R2"
            direction = 'rev' if read.is_reverse else 'fwd'
            label = f"{side}-{direction}"

            # Add it to the dict as a tuple of sequence and quality
            family_reads[
                families[read.query_name]
            ][
                label
            ].append(
                (read.query_alignment_sequence, read.query_alignment_qualities)
            )
            counter += 1

    print(f"Read in {counter:,} reads for {len(family_reads):,} families")
    return family_reads


def write_ssc(family_reads):
    """Write out the consensus for each part of each family, return a summary table."""

    # OUTPUT FILEPATHS
    FWD_R1 = f"{shard_ix}.FWD.R1.fastq.gz"
    REV_R1 = f"{shard_ix}.REV.R1.fastq.gz"
    FWD_R2 = f"{shard_ix}.FWD.R2.fastq.gz"
    REV_R2 = f"{shard_ix}.REV.R2.fastq.gz"

    # Make a list of the stats for each family
    ssc_summary = []

    # Number of families with data from only one strand
    n_single = 0

    # Open the output handles
    with gzip.open(FWD_R1, 'wt') as FWD_R1_HANDLE, \
        gzip.open(REV_R1, 'wt') as REV_R1_HANDLE, \
        gzip.open(FWD_R2, 'wt') as FWD_R2_HANDLE, \
        gzip.open(REV_R2, 'wt') as REV_R2_HANDLE:

        # Iteratively calculate the consensus sequences
        for cons_dict, stat in compute_ssc(family_reads):

            # If there were reads from both strands
            if cons_dict is not None:

                # Both of the -R1 and -R2 sequences must be the same length, and
                # the FWD and REV sequences cannot overhang the end of the molecules
                for d in ['fwd', 'rev']:
                    cons_dict[f"R1-{d}"], cons_dict[f"R2-{d}"] = trim_sscs(cons_dict[f"R1-{d}"], cons_dict[f"R2-{d}"])

                # Write out each read
                write_fastq(cons_dict["R1-fwd"], FWD_R1_HANDLE)
                write_fastq(cons_dict["R1-rev"], REV_R1_HANDLE)
                write_fastq(cons_dict["R2-fwd"], FWD_R2_HANDLE)
                write_fastq(cons_dict["R2-rev"], REV_R2_HANDLE)

                # Add the stats to the list
                ssc_summary.append(stat)

            # Otherwise
            else:

                # Increment the counter for single-stranded data
                n_single += 1

    print(f"SSCs with data from only one strand: {n_single:,}")

    # Convert the list of dicts to a DataFrame
    ssc_summary = pd.DataFrame(ssc_summary)

    print(f"Wrote out SSC data for {ssc_summary.shape[0]:,} families")

    # If no reads were written
    if ssc_summary.shape[0] == 0:

        # Delete the output files
        for fp in [FWD_R1, FWD_R2, REV_R1, REV_R2]:
            print(f"Removing output file {fp}")
            os.remove(fp)
            assert not os.path.exists(fp)

    # Return the family stats
    return ssc_summary


def compute_ssc(family_reads):
    """Compute the single-strand consensus for each set of reads."""

    # Iterate over each family
    for family_id, reads in family_reads.items():

        # Make a dict with the consensus sequences for each group
        group_consensus = dict()

        # Keep track of the stats for each group
        group_stats = dict(
            family=family_id
        )

        # Iterate over each subset of reads in the family
        for group_name, group_reads in reads.items():

            # Compute the consensus sequence
            cons_seq, cons_qual = compute_consensus(group_reads)

            # If the reads aligned in the reverse orientation
            if group_name.endswith("rev"):

                # Take the reverse complement of the sequence
                cons_seq = reverse_complement(cons_seq)

            # Add the consensus sequence to the dict
            group_consensus[group_name] = (family_id, cons_seq, cons_qual)

            # Add the stats for this group
            group_stats[f"{group_name}-n"] = len(group_reads)
            group_stats[f"{group_name}-len"] = len(cons_seq)

        # There must be 4 parts for each family
        if len(group_consensus) == 4:

            # Return the consensus sequences, as well as the stats
            yield group_consensus, group_stats

        else:

            # There must just be reads from one starnd
            assert len(group_consensus) == 2
            yield None, None


def compute_consensus(group_reads):
    """Compute the consensus nucleotide sequence and quality for a list of reads."""

    # If there is only one read
    if len(group_reads) == 1:

        # Then just return that read
        return group_reads[0][0], encode_quals(list(group_reads[0][1]))

    # Get the maximum read length
    max_rlen = max([len(r[0]) for r in group_reads])

    # Get the string representations of each read
    reads = [r[0] for r in group_reads]

    # For each read
    for i in range(len(reads)):

        # If the read is less than the maximum
        if len(reads[i]) < max_rlen:

            # Pad with N's
            reads[i] = reads[i] + "".join(["N" for _ in range(max_rlen - len(reads[i]))])

    # Format a multi-FASTA
    multi_fasta = "\\n".join(
        [
            f">{i}\\n{read}"
            for i, read in enumerate(reads)
        ]
    )

    # Convert to a Bio.Align object
    aln = SummaryInfo(
        AlignIO.read(
            StringIO(multi_fasta), 
            'fasta'
        )
    )

    # Compute the consensus sequence
    try:
        cons = aln.dumb_consensus(
            threshold=min_base_prop,
            ambiguous="N"
        )
    except Exception as e:
        print(multi_fasta)
        raise e

    # Remove any terminal N's
    cons = cons.rstrip("N")

    # Take the maximum quality score across all reads
    quals = defaultdict(list)
    for r in group_reads:
        for i, q in enumerate(r[1]):
            quals[i].append(q)
    quals = [max(quals[i]) for i in range(len(cons))]
    quals = encode_quals(quals)

    return str(cons), quals


def encode_quals(quals):
    return "".join([chr(n + 33) for n in quals])


def trim_sscs(r1, r2):
    """The reads from both strands must be the same length."""

    # Find the length of the shorter read
    min_len = min([len(r1[1]), len(r2[1])])

    # It's also very important that the length of the read does
    # not overhang the end of the molecule

    # First parse the size of the insert from the family name
    _, _, fwd_pos, rev_pos = r1[0].split("-")

    # Calculate the insert size
    insert_size = int(rev_pos) - int(fwd_pos)
    assert insert_size > 0, r1[0]

    # The read length cannot be longer than the insert size
    min_len = min([min_len, insert_size])

    # Trim down to the minimum length
    r1 = (r1[0], r1[1][:min_len], r1[2][:min_len])
    r2 = (r2[0], r2[1][:min_len], r2[2][:min_len])

    return r1, r2

def write_fastq(read_tuple, output_handle):
    """Write out a single read in FASTQ format."""

    family_id, seq, qual = read_tuple

    output_handle.write(
        f"@{family_id}\\n{seq}\\n+\\n{qual}\\n"
    )


def write_csv(family_stats, stats_csv):
    """Write out a list of dicts in CSV format."""

    # Get the header
    for k, v in family_stats.items():
        header = list(v.keys())
        break

    # Open the output file handle
    with gzip.open(stats_csv, "wt") as handle:

        # Write out the header
        handle.write(",".join(header) + '\\n')

        # Iterate over each item
        for r in family_stats:

            # Write out the row
            handle.write(",".join([
                str(r[cname])
                for cname in header
            ]) + '\\n')


# READ IN THE FAMILY ASSIGNMENT PER READ
families, family_size = read_family_assignments(families_csv_gz)

# PARSE THE READS BY FAMILY AND ORIENTATION
family_reads = parse_input_bam(bam, families)

# WRITE OUT THE SINGLE-STRAND CONSENSUS SEQUENCES
# ALSO GENERATE SUMMARY METRICS PER EACH SSC
ssc_summary = write_ssc(family_reads)

# IF DATA WAS FOUND
if ssc_summary.shape[0] > 0:

    # WRITE OUT THE STATS FOR EACH FAMILY
    ssc_summary.to_csv(stats_csv, index=None)
    print(f"Wrote out to {stats_csv}")

print("DONE")