#!/usr/bin/env python3

import argparse
import gzip
import logging
from Bio.Seq import reverse_complement
import pysam
import pandas as pd


class TrimOverhang:
    """Useful object to trim the overhanging 3' bases from paired-end reads."""

    def __init__(
        self,
        input_bam=None,
        input_positions=None,
        output_read1=None,
        output_read2=None
    ):

        # Set up logging
        self.init_logging()

        # Add the kwargs to the object namespace
        self.input_bam = input_bam
        self.logger.info(f"Input BAM: {self.input_bam}")
        self.input_positions = input_positions
        self.logger.info(f"Input positions: {self.input_positions}")
        self.output_read1 = output_read1
        self.logger.info(f"Output FASTQ (R1): {self.output_read1}")
        self.output_read2 = output_read2
        self.logger.info(f"Output FASTQ (R2): {self.output_read2}")

        # Get the lengths of every sequenced molecule
        self.get_insert_length()

        # Stream through the BAM and write out trimmed reads
        self.trim_reads()

    def init_logging(self):
        """Set up logging."""

        # Set the level of the logger to INFO
        logFormatter = logging.Formatter(
            '%(asctime)s %(levelname)-8s [TrimOverhang] %(message)s'
        )
        self.logger = logging.getLogger('TrimOverhang')
        self.logger.setLevel(logging.INFO)

        # Write to STDOUT
        consoleHandler = logging.StreamHandler()
        consoleHandler.setFormatter(logFormatter)
        self.logger.addHandler(consoleHandler)

    def get_insert_length(self):
        """Get the lengths of every sequenced molecule."""
        
        # Read in the table listing the location of each pair of reads
        positions = pd.read_csv(self.input_positions)

        self.logger.info(f"Read in {self.input_positions} - {positions.shape[0]:,} individual reads.")
        self.logger.info(positions.head())

        # Pivot the table to just contain the outermost coordinates of each read pair
        positions = positions.pivot_table(
            index='id',
            columns='direction',
            values='pos'
        )

        self.logger.info(f"Pivoted to yield outermost alignment coordinates - {positions.shape[0]:,} read pairs")
        self.logger.info(positions.head())

        # Drop any data where only 1/2 reads is aligned
        positions = positions.dropna().applymap(int)
        self.logger.info(positions.head())

        self.logger.info(f"Removed any unpaired reads - {positions.shape[0]:,} read pairs remaining")

        # Compute the length of the sequenced molecule
        positions = positions.assign(
            total_length = positions.apply(
                lambda r: r['rev'] - r['fwd'] + 1,
                axis=1
            )
        )

        self.logger.info("Computed fragment lengths")
        self.logger.info(positions.head())

        # Remove any pairs with implausible fragment lengths
        positions = positions.query(
            "total_length > 1"
        ).query(
            "total_length < 10000"
        )

        self.logger.info(f"Filtered out fragments <=1 and >= 10000: {positions.shape[0]:,} read pairs remaining")
        self.logger.info(positions.head())

        # Make a dictionary with the read length values and attach it to the object
        self.lengths = positions['total_length'].to_dict()

    def trim_reads(self):
        """Stream through the BAM and write out trimmed reads."""

        # Keep a buffer of reads so that we can write each pair
        # of reads in sync across both files
        self.read_buffer = dict(
            R1=dict(),
            R2=dict()
        )

        # Open the BAM file for reading
        self.logger.info(f"Opening {args.input_bam} for reading")
        bam = pysam.AlignmentFile(args.input_bam, "rb")

        # Open the FASTQ files for writing
        self.logger.info(f"Opening {args.output_read1} for writing")
        self.R1 = gzip.open(args.output_read1, "wt")
        self.logger.info(f"Opening {args.output_read2} for writing")
        self.R2 = gzip.open(args.output_read2, "wt")

        # Keep track of the read pairs which were written
        self.reads_written = set([])

        # Iterate over every read in the alignment
        for read in bam:

            # Write out a read if its mate is already found,
            # otherwise add it to the buffer
            self.process_read(read)

        # Close the file handles
        self.logger.info("Closing all file handles")
        bam.close()
        self.R1.close()
        self.R2.close()

        # Report how many read pairs were written
        self.logger.info(f"Wrote out {len(self.reads_written):,} read pairs")

        # Report how many reads are left in the buffer
        for read_direction, read_set in self.read_buffer.items():

            self.logger.info(f"Number of unpaired {read_direction} reads which were not written: {len(read_set):,}")

        # Make sure that the expected set of reads were written
        assert len(self.reads_written) >= len(self.lengths) * 0.9, f"Did not write out the full set of {len(self.lengths):,} read pairs"

    def process_read(self, read):
        """Write out a read if its mate is already found, otherwise add it to the buffer."""

        # Get the length of the fragment, inferred from the outermost aligned positions
        read_len = self.lengths.get(read.query_name)

        # If this read was not part of a proper pair
        if read_len is None:

            # Skip it
            return

        # If this read was already written
        if read.query_name in self.reads_written:

            # Skip it
            return

        # Get the read pair from the buffer, if any exists
        read_mate = self.get_mate(read)

        # If the mate is in the buffer
        if read_mate is not None:

            # Write out both reads
            self.write_read_fastq(read)
            self.write_read_fastq(read_mate)

            # Remove the mate read from the buffer
            self.clear_mate(read)

            # Mark this read as having been written out as a pair
            self.reads_written.add(read.query_name)

        # If the mate is not in the buffer
        else:

            # Add the read to the buffer
            self.add_to_buffer(read)

    def add_to_buffer(self, read):
        """Add a read to the buffer."""

        self.read_buffer[
            'R1' if read.is_read1 else 'R2'
        ][
            read.query_name
        ] = read

    def get_mate(self, read):
        """Return the mate of a read from the buffer, if any."""

        return self.read_buffer[
            'R2' if read.is_read1 else 'R1'
        ].get(
            read.query_name
        )
        
    def clear_mate(self, read):
        """Remove the mate read from the buffer."""

        del self.read_buffer[
            'R2' if read.is_read1 else 'R1'
        ][
            read.query_name
        ]
        
    def write_read_fastq(self, read):
        """Write out a single read in FASTQ format, trimming if it overhangs the end."""

        # Get the sequence of the read and its quality
        read_name = read.query_name
        read_nucl = read.query_alignment_sequence
        read_qual = self.encode_quals(read.query_alignment_qualities)

        # Get the barcode assigned to the read
        barcode = read.get_tag("BC")

        # Get the length of the molecule, inferred from the paired read
        read_len = self.lengths[read_name]

        # Both lengths must match
        assert len(read_nucl) == len(read_qual)

        # If the read is longer than the fragment
        if len(read_nucl) > read_len:

            # Then trim both strings
            read_nucl = read_nucl[:read_len]
            read_qual = read_qual[:read_len]

        # Set the output handle depending on whether the read is R1 or R2
        output_handle = self.R1 if read.is_read1 else self.R2

        # If the read is mapped on the reverse strand
        if read.is_reverse:

            # Reverse complement the nucleotide sequence
            read_nucl = reverse_complement(read_nucl)

            # Reverse the quality
            read_qual = read_qual[::-1]

        # Write out in FASTQ format
        output_handle.write(
            f"@{read_name} BC:Z:{barcode}\n{read_nucl}\n+\n{read_qual}\n"
        )

    def encode_quals(self, quals):
        return "".join([chr(n + 33) for n in quals])

# If this file is being executed as a script
if __name__ == "__main__":

    # Parse the user-provided arguments
    parser = argparse.ArgumentParser(
        description="Trim the overhang from any aligned reads"
    )

    parser.add_argument(
        '--input-bam',
        help='Input BAM file containing reads to be trimmed.',
        required=True
    )
    parser.add_argument(
        '--input-positions',
        help='Input CSV file listing the positions of each read.',
        required=True
    )
    parser.add_argument(
        '--output-read1',
        help='Output reads (R1) in FASTQ format.',
        required=True
    )
    parser.add_argument(
        '--output-read2',
        help='Output reads (R2) in FASTQ format.',
        required=True
    )
    args = parser.parse_args()

    # Trim the reads and write to the output
    # Map the command line arguments directly to the kwargs
    TrimOverhang(
        **args.__dict__
    )
