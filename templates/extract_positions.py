#!/usr/bin/env python3

import gzip
import os
import pysam

bam_fp = "${bam}"
print(f"Reading input from {bam_fp}")
assert os.path.exists(bam_fp)
output_fp = "read_positions.csv.gz"
print(f"Writing output to {output_fp}")

# Define the columns in the CSV output
header = ['id', 'direction', 'strand', 'chr', 'pos', 'barcode']
def write_csv(d, handle):
    """Write out one line of the CSV."""
    handle.write(
        ",".join([
            str(d[cname])
            for cname in header
        ]) + "\\n"
    )

counter = 0

# Open the input and output files
with pysam.AlignmentFile(bam_fp, "rb") as bam, gzip.open(output_fp, "wt") as csv:

    # Write the CSV header
    csv.write(",".join(header) + "\\n")

    # Iterate over each read
    for read in bam:

        # Write out the information for this read
        write_csv(
            dict(
                # The ID of the read
                id=read.query_name,
                # strand = R1 / R2
                strand = "R1" if read.is_read1 else "R2",
                # direction = fwd / rev
                direction = "rev" if read.is_reverse else "fwd",
                # The position on the contig (transform to 1-based indexing)
                pos = read.reference_end + 1 if read.is_reverse else read.reference_start + 1,
                # The name of the contig that the read is aligned to
                chr = read.reference_name,
                # The barcode sequence
                barcode = read.get_tag("BC")
            ),
            csv
        )
        counter += 1

print(f"Wrote out information for {counter:,} reads")

