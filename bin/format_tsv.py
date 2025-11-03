#!/usr/bin/env python3

from collections import defaultdict
import gzip
import sys

pileup_in = sys.argv[1]
tsv_out = sys.argv[2]

print(f"Reading in {pileup_in}")
print(f"Writing out {tsv_out}")

# Define the columns of the output
output_cols = [
    "chrom",
    "pos",
    "ref",
    "depth",
    "mutation_count",
    "T",
    "C",
    "G",
    "A",
    "mutation_classes"
]

def parse_line(line):
    """Parse a single line from a pileup file."""

    fields = dict(zip(
        [
            "chrom",
            "pos",
            "ref",
            "depth",
            "bases",
            "qual",
        ],
        line.rstrip("\n").split("\t")
    ))

    # Parse the number of mutations from the 'base' string
    for k, v in parse_muts(fields['bases']).items():
        fields[k] = v

    return "\t".join([
        str(fields[cname])
        for cname in output_cols
    ]) + "\n"


def parse_muts(base_string):
    """Parse the number of mutations from the base string."""

    muts = dict(
        A=0,
        T=0,
        C=0,
        G=0
    )

    mut_counts = 0

    for base in base_string:
        if muts.get(base) is not None:
            muts[base] += 1
            mut_counts += 1

    muts["mutation_classes"] = sum([v > 0 for v in muts.values()])
    muts["mutation_count"] = mut_counts

    return muts

with gzip.open(pileup_in, "rt") as handle_in, gzip.open(tsv_out, "wt") as handle_out:

    # Write the header
    handle_out.write("\t".join(output_cols) + "\n")

    for line in handle_in:
        handle_out.write(
            parse_line(
                line
            )
        )