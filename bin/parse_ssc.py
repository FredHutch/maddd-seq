#!/usr/bin/env python3
"""
Parse the SSC data in order to:
  - Call adducts from the SSC data
  - Call SNPs from the SSC data
  - Count the number of mutations and adducts per DSC
  - Make a table of all mutations and adducts per DSC
  - Write out a DSC BAM and SSC BAM for each possible
    level of filtering based on the maximum number of
    allowable mutations and adducts per DSC
  - Make a summary table of the total number of adducts
    and SNPs for each of the possible levels of filtering
"""

from collections import defaultdict
import gzip
import json
import logging
import os
import pandas as pd
import pysam
import sys

# Set up logging
logFormatter = logging.Formatter(
    '%(asctime)s %(levelname)-8s [parse_ssc] %(message)s'
)
logger = logging.getLogger()
logger.setLevel(logging.INFO)

# Write to STDOUT
consoleHandler = logging.StreamHandler()
consoleHandler.setFormatter(logFormatter)
logger.addHandler(consoleHandler)

# Parse the input arguments
specimen = sys.argv[1]
logger.info(f"Processing specimen: {specimen}")

# Input filepaths for filtered SSC sequences
input_pos_bam = "POS.SSC.bam"
assert os.path.exists(input_pos_bam)
input_neg_bam = "NEG.SSC.bam"
assert os.path.exists(input_neg_bam)


def merge_strands(pos_ssc, neg_ssc):
    """Merge information from both strands to identify SNPs and adducts."""

    # Output data is a dictionary keyed by family
    output = dict()

    # Also keep track of the number of mutations and adducts
    # as a function of the position in the read
    base_positions = dict(
        nreads=defaultdict(int),
        variants=defaultdict(int),
        adducts=defaultdict(int)
    )

    # Keep a list of adducts
    adducts = []

    # Keep a list of all families which contain adducts
    adduct_families = []

    # Iterate over each family
    for family_id in pos_ssc:

        # Make sure that the family is in the negative strand data as well
        assert neg_ssc.get(family_id) is not None, f"Missing negative strand information for {family_id}"

        # Merge data from fwd/rev, pos/neg
        output[family_id] = merge_single_family(pos_ssc.get(family_id), neg_ssc.get(family_id))

        # If there are any adducts
        if len(output[family_id]['adducts']) > 0:

            # Add this family to the list
            adduct_families.append(family_id)

            # Add any adducts to the list
            add_adduct_info(
                output[family_id],
                adducts
            )

        # Add the nreads/snps/adducts data as a function of read position
        add_base_position_info(
            output[family_id],
            pos_ssc.get(family_id),
            neg_ssc.get(family_id),
            base_positions
        )

        # Add the number of aligned bases
        output[family_id]['nbases'] = calc_nbases(pos_ssc.get(family_id), neg_ssc.get(family_id))

    return output, base_positions, adducts, adduct_families


def get_read_pos(ref_pos, pos_family):
    """Return the position of a reference base, relative to the read."""
    return min(
        abs(pos_family['fwd']['start'] - ref_pos),
        abs(pos_family['rev']['start'] - ref_pos)
    )

def add_base_position_info(merged_family, pos_family, neg_family, base_positions):
    """Add the nreads/snps/adducts data as a function of read position."""

    # Increment over all of the positions for which adducts were detected
    for adduct_pos in merged_family['adducts']:
        # Get the position of the mutation inside the read
        read_pos = get_read_pos(adduct_pos, pos_family)
        # Increment the counter
        base_positions['adducts'][read_pos] += 1

    # Increment over all of the positions for which mutations were detected
    for adduct_pos in merged_family['mutations']:
        # Get the position of the mutation inside the read
        read_pos = get_read_pos(adduct_pos, pos_family)
        # Increment the counter
        base_positions['snps'][read_pos] += 1

    # Iterate over both the forward and reverse reads
    for read_dir in ['fwd', 'rev']:

        # Calculate the read length
        rlen = abs(pos_family[read_dir]['start'] - pos_family[read_dir]['end'])

        # Iterate over the read length
        for base_pos in range(rlen):

            # Increment the counter
            base_positions['nreads'][base_pos] += 1


def calc_nbases(pos_family, neg_family):
    """
    Calculate the number of aligned bases.
    The positive and negative strands have been trimmed to the same
    length, and so we only have to do the calculation once.
    """

    # If the reads overlap
    if pos_family['fwd']['end'] >= pos_family['rev']['start']:

        # Calculate the total span
        return pos_family['rev']['end'] - pos_family['fwd']['start']

    # If the reads do not overlap
    else:

        # Sum the span from each end
        return (pos_family['fwd']['end'] - pos_family['fwd']['start']) + (pos_family['rev']['end'] - pos_family['rev']['start'])


def merge_single_family(pos_family, neg_family, allowed_bases=set(['A', 'T', 'C', 'G'])):
    """Merge the information from a single family."""

    # Keep track of mutations and adducts
    output = dict(
        ref_name=pos_family['fwd']['ref_name'],
        mutations=dict(),
        adducts=dict(),
        references=dict()
    )

    # Get a list of all positions which have variants
    var_pos_list = list(
        set(pos_family['fwd']['variants'].keys()) | \
        set(pos_family['rev']['variants'].keys()) | \
        set(neg_family['fwd']['variants'].keys()) | \
        set(neg_family['rev']['variants'].keys())
    )

    # Get the corresponding reference base at each of those positions
    reference_bases = dict()
    for family_dict in [pos_family, neg_family]:
        for dir_dict in family_dict.values():
            for pos, base in dir_dict.get('references', {}).items():
                reference_bases[pos] = base.upper()

    # Iterate over each position
    for var_pos in var_pos_list:

        # For each of the parts of this family, classify this position
        family_dat = dict(
            fwd_pos=classify_position(var_pos, pos_family['fwd']),
            rev_pos=classify_position(var_pos, pos_family['rev']),
            fwd_neg=classify_position(var_pos, neg_family['fwd']),
            rev_neg=classify_position(var_pos, neg_family['rev'])
        )

        # If we are lacking base calls for any aligned regions
        if not all(map(lambda b: b in ['ref', 'unaligned'] or b in allowed_bases, family_dat.values())):

            # Skip the position
            continue

        # If the fwd/pos read is aligned
        if family_dat['fwd_pos'] != 'unaligned':
            
            # And the rev/pos read is not consistent
            if family_dat['rev_pos'] not in ['unaligned', family_dat['fwd_pos']]:

                # Skip it
                continue

            # Otherwise, set the positive strand base from the forward read
            pos_base = family_dat['fwd_pos']

        # If the fwd/pos base is not aligned
        else:

            # Set the positive strand base from the reverse read
            pos_base = family_dat['rev_pos']

        # If the fwd/neg read is aligned
        if family_dat['fwd_neg'] != 'unaligned':
            
            # And the rev/neg read is not consistent
            if family_dat['rev_neg'] not in ['unaligned', family_dat['fwd_pos']]:

                # Skip it
                continue

            # Otherwise, set the negative strand base from the forward read
            neg_base = family_dat['fwd_neg']

        # If the fwd/neg base is not aligned
        else:

            # Set the negative strand base from the reverse read
            neg_base = family_dat['rev_neg']

        # If either base is unaligned
        if pos_base == "unaligned" or neg_base == "unaligned":

            # Skip the position
            continue

        # If both strands agree
        if pos_base == neg_base:

            # Then this is a mutation
            output['mutations'][var_pos] = pos_base

        # if the strands do not agree, and the positive strand is unchanged
        elif pos_base == 'ref':
            # Add the adduct to the list
            output['adducts'][var_pos] = neg_base, 'neg'

        # If the negative strand is unchanged
        elif neg_base == 'ref':
            # Add the adduct to the list
            output['adducts'][var_pos] = pos_base, 'pos'

        # If neither strand matches the reference
        else:

            # The positive strand is a SNP
            output['mutations'][var_pos] = pos_base

            # The negative strand is an adduct
            output['adducts'][var_pos] = neg_base, 'neg'

        # Keep track of the reference base
        output['references'][var_pos] = reference_bases[var_pos]

    return output

def classify_position(var_pos, family_info):
    """For a single position, characterize this read as being 'ref', 'unaligned', or the variant sequence."""

    # If this position is outside the aligned region
    if var_pos < family_info['start'] or var_pos >= family_info['end']:
        return 'unaligned'

    # Otherwise, return the variant base, substituting
    # with 'ref' if the position does not have a variant noted
    else:
        return family_info['variants'].get(var_pos, 'ref')


def complement(base):
    """Return the complementary base."""

    return dict(
        A='T',
        T='A',
        C='G',
        G='C'
    )[base]


def iupac(base1, base2):
    """https://www.bioinformatics.org/sms/iupac.html"""

    return dict(
        A=dict(
            A="A",
            T="W",
            C="M",
            G="R",
        ),
        T=dict(
            A="W",
            T="T",
            C="Y",
            G="K",
        ),
        C=dict(
            A="M",
            T="Y",
            C="C",
            G="S",
        ),
        G=dict(
            A="R",
            T="K",
            C="S",
            G="G",
        ),
    ).get(base1, {}).get(base2, "N")


class ParseSSC:
    """Class used to analyze SSC data from BAM inputs."""

    def __init__(self, specimen):

        # Record the specimen name
        self.specimen = specimen

        # Keep track of the reference sequence
        # Key by ref_name and position
        self.refseq = defaultdict(lambda: dict())

        # Record information from each read
        self.read_info = defaultdict(
            # Keep track of data from the positive and negative strands,
            # for the forward and reverse reads
            lambda: {
                strand: {
                    orient: dict()
                    for orient in ["fwd", "rev"]
                }
                for strand in ["pos", "neg"]
            }
        )

        # Keep information about the position of variants
        # relative to the read start
        self.base_positions = {
            key: defaultdict(int)
            for key in ['adducts', 'variants', 'nreads']
        }

        # Parse the information from the positive strand
        self.parse_bam(fp=input_pos_bam, strand="pos")
        # Parse the information from the negative strand
        self.parse_bam(fp=input_neg_bam, strand="neg")

        # Merge the information from the forward and reverse
        # reads of each SSC
        self.merge_fwd_rev_per_strand()

        # Merge the information from both strands
        self.merge_pos_neg_strands()

        # Write out the total information without filtering
        self.write_output(folder="all", max_vars=None)

        # For each of the unique values of the total number of variants and adducts per read
        for max_vars in list(set([
            self.dsc_info[family_id]["total_variants_and_adducts"]
            for family_id in self.dsc_info
        ])):

            # Write out the total information, filtering to that maximum number
            self.write_output(folder=f"max_variants_{max_vars}", max_vars=max_vars)

    def parse_bam(self, fp=None, strand=None):
        """Parse all of the data from a single BAM."""

        assert strand in ["pos", "neg"], f"Unrecognized strand '{strand}'"

        # Open the input BAM file for reading
        logger.info(f"Reading from {fp}")
        with pysam.AlignmentFile(fp, "rb") as bam:

            # For each of the reads on the positive strand    
            for i, read in enumerate(bam):

                if i >= 10000:
                    break

                # Get the stats for this read
                family_id, orientation, read_stats = self.parse_read(read)

                # Add to the dataset
                self.read_info[
                    family_id
                ][
                    strand
                ][
                    orientation
                ] = read_stats

    def parse_read(self, read):
        """
        For a single read, capture the following information:
        - start: leftmost position of alignment
        - end: rightmost position of alignment
        - variants: dict with position/base for any base that differs from the reference
        """

        # Family ID
        family_id = read.query_name
        
        # Orientation
        orient = 'rev' if read.is_reverse else 'fwd'
        
        # Details for the position and variants in the read
        read_details = dict(
            # Chromosome / contig name
            ref_name = read.reference_name,
            # Leftmost position
            start = read.reference_start,
            # Rightmost position
            end = read.reference_end,
            # Variants
            variants = self.parse_variants(read)
        )

        # Tuple, Family ID, orientation, and a dict with position and variants
        return family_id, orient, read_details

    def parse_variants(self, read, allowed_nucs=set(['A', 'T', 'C', 'G'])):
        """
        For any position in which the read differs from the reference,
        include the reference position and the variant base in a dict.
        Skip any position which does not contain one of the `allowed_nucs`.
        Lowercase bases in the reference will automatically be transformed
        into uppercase. Based on this behavior, soft-masked bases will be
        included in all of the mutational positions.
        """

        # Encode the as a dict of positions where the read does not match the reference
        variants = dict()

        # Get the reference sequence
        rseq = read.get_reference_sequence()

        # get_aligned_pairs() returns a tuple of positions
        for qpos, rpos in read.get_aligned_pairs():

            # If there is an indel
            if rpos is None or qpos is None:
                
                # Skip it
                continue

            rpos = int(rpos)
            qpos = int(qpos)

            # Get the aligned and reference bases
            qbase = read.query_sequence[qpos].upper()
            rbase = rseq[rpos - read.reference_start].upper()

            # Record the reference position
            self.refseq[read.reference_name][rpos] = rbase

            # If the reference base has been masked
            if rbase not in allowed_nucs or qbase not in allowed_nucs:

                # Skip it
                continue

            # If the query matches the reference
            if qbase == rbase:

                # Skip it
                continue

            # Otherwise, record the query base in the output
            variants[rpos] = qbase

        return variants

    def merge_fwd_rev_per_strand(self):
        """Merge the information for each forward and reverse read per strand."""

        # Record the information for each strand
        # Keyed first by family ID
        self.ssc_info = defaultdict(
            lambda: {
                # Keyed second by pos/neg
                strand: dict(
                    # Chromosome / contig name
                    ref_name = None,
                    # Leftmost position
                    start = None,
                    # Rightmost position
                    end = None,
                    # Consensus sequence
                    cons = None,
                    # Variants relative to the reference
                    variants = dict()
                    # Primary key is reference position, values are:
                    # dict(
                    #     readpos = (position relative to read)
                    #     var     = (variant base)
                    #     ref     = (reference base)
                    # )
                )
                for strand in ["pos", "neg"]
            }
        )

        # For each family
        for family_id, family_reads in self.read_info.items():

            # For each strand
            for strand, strand_reads in family_reads.items():

                # Merge a single pair of forward and reverse reads
                self.merge_read_pair(family_id, strand, strand_reads)

    def merge_read_pair(self, family_id, strand, strand_reads):
        """Merge the forward and reverse reads for a single strand of a single family."""

        # If we don't have forward and reverse reads
        if "fwd" not in strand_reads or "rev" not in strand_reads:
            logger.info("Unexpected - didn't find forward and reverse")
            logger.info(json.dumps(strand_reads))
            return

        # If the forward and reverse strands are on different chromosomes
        ref_name = strand_reads["fwd"]["ref_name"]
        if ref_name != strand_reads["rev"].get("ref_name"):

            # Log and stop
            logger.info(f"Unexpected -- reads are on different references ({family_id} - {strand})")
            logger.info(json.dumps(strand_reads))
            return

        # If the forward read is not 5' to the reverse read
        if strand_reads["fwd"]["start"] >= strand_reads["rev"]["end"]:

            # Log and stop
            logger.info(f"Unexpected -- reads are not oriented inwards ({family_id} - {strand})")
            logger.info(json.dumps(strand_reads))
            return

        # Record the positional information
        self.ssc_info[family_id][strand]["ref_name"] = ref_name
        self.ssc_info[family_id][strand]["start"] = strand_reads["fwd"]["start"]
        self.ssc_info[family_id][strand]["end"] = strand_reads["rev"]["end"]

        # Get the set of positions covered by the read pair
        covered_bases = dict()

        # Combine the variants from the forward and reverse reads
        for fwd_rev, strand_info in strand_reads.items():

            # Add the bases covered by this read
            covered_bases[fwd_rev] = set(range(strand_info["start"], strand_info["end"]+1))

            # Iterate over the variant bases
            for refpos, variant_base in strand_info["variants"].items():

                # Get the position relative to the read
                if fwd_rev == "fwd":
                    readpos = (refpos - strand_info["start"]) + 1
                else:
                    readpos = (strand_info["end"] - refpos) + 1

                # Add the data to the collection of variants
                self.ssc_info[family_id][strand]["variants"][refpos] = dict(
                    readpos=readpos,
                    var=variant_base,
                    ref=self.refseq[ref_name][refpos]
                )

        # Get the combined set of covered positions
        covered_bases = covered_bases["fwd"] | covered_bases["rev"]

        # Build the consensus sequence
        cons = []

        # Iterate over each position
        for pos in range(
            self.ssc_info[family_id][strand]["start"],
            self.ssc_info[family_id][strand]["end"] + 1
        ):

            # If the position is covered
            if pos in covered_bases:

                # If the position is a variant, add it
                # otherwise add the reference
                cons.append(
                    self.ssc_info[family_id][strand]["variants"].get(
                        pos,
                        {}
                    ).get(
                        "var",
                        self.refseq[ref_name].get(pos, "N")
                    )
                )

            # If the position is not covered
            else:

                # Add N
                cons.append("N")

        # Concatenate and save the consensus sequence
        self.ssc_info[family_id][strand]["cons"] = "".join(cons)

    def merge_pos_neg_strands(self):
        """Merge the information for the positive and negative strands per family."""

        # Record the information for each strand
        # Keyed by family ID
        self.dsc_info = defaultdict(
            lambda: dict(
                # Chromosome / contig name
                ref_name = None,
                # Leftmost position
                start = None,
                # Rightmost position
                end = None,
                # Consensus sequence
                cons = None,
                # Number of bases sequenced
                nbases = 0,
                # Adducts (single strand)
                # {
                #   pos: dict(
                #     strand=pos/neg,
                #     var=base,  <- from the appropriate strand
                #     ref=base,  <- from the appropriate strand
                #   )
                # }
                adducts = dict(),
                # Variants (double stranded)
                # {
                #   pos: dict(
                #     var=base,  <- from the positive strand
                #     ref=base,  <- from the positive strand
                #   )
                # }
                variants = dict(),
                # Keep track of the total number of variants and adducts
                total_variants_and_adducts = 0
            )
        )
        
        # For each family
        for family_id, family_strands in self.ssc_info.items():

            # Merge the strands
            self.merge_strands(family_id, family_strands)

    def merge_strands(self, family_id, family_strands, allowed_nucs=set(['A', 'T', 'C', 'G'])):
        """Merge the information for both strands."""

        # Get the inner positions for the start and stop
        start_pos = max(family_strands["pos"]["start"], family_strands["neg"]["start"])
        end_pos = min(family_strands["pos"]["end"], family_strands["neg"]["end"])

        self.dsc_info[family_id]["start"] = start_pos
        self.dsc_info[family_id]["end"] = end_pos

        # Trim the consensus sequences for each strand
        strand_cons = {
            strand: strand_dict["cons"][
                strand_dict["start"] - start_pos: len(strand_dict["cons"]) - (strand_dict["end"] - end_pos)
            ].upper()
            for strand, strand_dict in family_strands.items()
        }

        # Build a list of the double-stranded consensus sequence
        cons_list = []

        # Get the reference name
        ref_name = family_strands["pos"]["ref_name"]
        self.dsc_info[family_id]["ref_name"] = ref_name

        # Iterate over the sequence at each position
        for refpos, pos_base, neg_base in zip(range(start_pos, end_pos + 1), strand_cons["pos"], strand_cons["neg"]):

            # If either strand was not sequenced
            if pos_base not in allowed_nucs or neg_base not in allowed_nucs:

                # Add an N to the consensus
                cons_list.append("N")

                # Move on
                continue

            # Get the reference base at this position
            refbase = self.refseq[ref_name][refpos]

            # Get the shortest distance to the either end
            readpos = min(refpos - start_pos, end_pos - refpos) + 1
            assert readpos > 0, (refpos, start_pos, end_pos)

            # Increment the counter with the number of bases sequenced
            self.dsc_info[family_id]["nbases"] += 1
            self.base_positions["nreads"][readpos] += 1

            # Add the merged base to the consensus
            cons_list.append(iupac(pos_base, neg_base))

            # If the reference is not known at this position
            if refbase not in allowed_nucs:

                # Move on
                continue

            # If both bases are mismatched
            if pos_base != refbase and neg_base != refbase:

                # If they are the same
                if pos_base == neg_base:

                    # It is a variant
                    self.dsc_info[family_id]["variants"][refpos] = dict(
                        var=pos_base,
                        ref=refbase
                    )
                    self.base_positions["variants"][readpos] += 1

                # If they are different
                else:

                    # The positive strand is the variant
                    self.dsc_info[family_id]["variants"][refpos] = dict(
                        var=pos_base,
                        ref=refbase
                    )
                    self.base_positions["variants"][readpos] += 1

                    # And the negative strand is the adduct
                    self.dsc_info[family_id]["adducts"][refpos] = dict(
                        strand="neg",
                        var=complement(neg_base),
                        ref=complement(refbase)
                    )
                    self.base_positions["adducts"][readpos] += 1

                # Increment the total number of variants and adducts
                self.dsc_info[family_id]["total_variants_and_adducts"] += 1

            # If only the positive strand is mismatched
            elif pos_base != refbase:

                # The positive strand is the adduct
                self.dsc_info[family_id]["adducts"][refpos] = dict(
                    strand="pos",
                    var=pos_base,
                    ref=refbase
                )
                self.base_positions["adducts"][readpos] += 1

                # Increment the total number of variants and adducts
                self.dsc_info[family_id]["total_variants_and_adducts"] += 1

            # If only the negative strand is mismatched
            elif neg_base != refbase:

                # The negative strand is the adduct
                self.dsc_info[family_id]["adducts"][refpos] = dict(
                    strand="neg",
                    var=complement(neg_base),
                    ref=complement(refbase)
                )
                self.base_positions["adducts"][readpos] += 1

                # Increment the total number of variants and adducts
                self.dsc_info[family_id]["total_variants_and_adducts"] += 1

        self.dsc_info[family_id]["cons"] = "".join(cons_list)

    def write_output(self, folder=None, max_vars=None):
        """Write all outputs to a folder, optionally filtering by total number of variants."""

        # If the folder doesn't exist
        if not os.path.exists(folder):

            # Create it
            logger.info(f"Creating folder {folder}")
            os.mkdir(folder)

        # Make a set of all families which do not exceed the filter
        keep_families = set([
            family_id
            for family_id, dsc in self.dsc_info.items()
            if max_vars is None or dsc["total_variants_and_adducts"] <= max_vars
        ])

        # Save the adduct information as GTF
        self.write_adducts_gtf(folder=folder, keep_families=keep_families)

        # Save the list of all families which contain adducts
        self.write_adduct_family_list(folder=folder, keep_families=keep_families)

        # Save the total information to JSON
        self.write_total_json(folder=folder, keep_families=keep_families)

        # Format all of the output, both by contig and overall
        summary_dat, by_chr, snp_base_changes, adduct_base_changes = self.format_summary(keep_families=keep_families)

        # Save the summary information to JSON
        summary_json_fpo = os.path.join(folder, "summary.json.gz")
        logger.info(f"Writing summary information to {summary_json_fpo}")
        with gzip.open(summary_json_fpo, "wt") as handle:
            json.dump(summary_dat, handle)

        # Output path for the data grouped by contig
        by_chr_fpo = os.path.join(folder, f"{specimen}.by_chr.csv.gz")
        logger.info(f"Writing data grouped by contig to {by_chr_fpo}")

        # Save the by-chr information to CSV
        by_chr.to_csv(by_chr_fpo, index_label='chr')

        # Output path for SNP data grouped by base change
        snp_base_fpo = os.path.join(folder, f"{specimen}.snps_by_base.csv.gz")
        logger.info(f"Writing SNP data grouped by base change to {snp_base_fpo}")

        # Save the SNP by-base information to CSV
        snp_base_changes.to_csv(snp_base_fpo, index_label='base')

        # Output path for adduct data grouped by base change
        adduct_base_output = os.path.join(folder, f"{specimen}.adducts_by_base.csv.gz")
        logger.info(f"Writing adduct data grouped by base change to {adduct_base_output}")

        # Save the adduct by-base information to CSV
        adduct_base_changes.to_csv(adduct_base_output, index_label='base')

        # Output path for variant data organized by position in the reads
        base_positions_output = os.path.join(folder, f"{specimen}.by_read_position.csv.gz")
        logger.info(f"Output path: {base_positions_output}")

        # Format the data by read position as a DataFrame
        base_positions = pd.DataFrame(self.base_positions).fillna(0).applymap(int)
        base_positions.to_csv(base_positions_output, index_label="pos")

        # Write out the filtered DSC and SSC as BAM
        for info, prefix, flag in [
            (self.dsc_info, "DSC", 99),
            # ({family_id: ssc["pos"] for family_id, ssc in self.ssc_info.items()}, "SSC.POS", 99),
            # ({family_id: ssc["neg"] for family_id, ssc in self.ssc_info.items()}, "SSC.NEG", 83)
        ]:
            self.write_bam(
                info,
                fp=os.path.join(folder, f"{prefix}.max_variants_{max_vars}.bam"),
                keep_families=keep_families,
                flag=flag
            )

    def write_bam(self, info, keep_families=None, fp=None, flag=99):
        """
        Write out a BAM file with all of the reads from `info`
        `info` is a dict with 
        """

        logger.info(f"Copying header from {input_pos_bam}")
        with pysam.AlignmentFile(input_pos_bam, "r") as template:

            # Map each reference name to an id
            reference_id_map = {
                d["SN"]: i
                for i, d in enumerate(template.header.to_dict()["SQ"])
            }

            logger.info(f"Writing out BAM to {fp}")
            with pysam.AlignmentFile(fp, "wb", template=template) as outf:

                # Iterate over each family
                for family_id, family_dat in info.items():
                    
                    # If the family is in the list to keep
                    if family_id in keep_families and family_dat["ref_name"] in reference_id_map:

                        seqlen = (family_dat["end"] - family_dat["start"]) + 1
                        a = pysam.AlignedSegment()
                        a.query_name = family_id
                        a.query_sequence = family_dat["cons"]
                        a.flag = flag
                        a.reference_id = reference_id_map[family_dat["ref_name"]]
                        a.reference_start = family_dat["start"]
                        a.mapping_quality = 20
                        a.cigar = [(0,seqlen)]
                        a.query_qualities = pysam.qualitystring_to_array("".join(["?" for _ in family_dat["cons"]]))
                        outf.write(a)

    def format_summary(self, keep_families=None):
        """Summarize the output, both by contig and overall."""

        # Count up the number of molecules, bases, snps, and adducts
        # both overall
        total_counts = defaultdict(int)
        # by chromosome
        chr_counts = defaultdict(lambda: defaultdict(int))
        # by the base change for mutations
        # A -> T, T -> C, etc.
        snp_base_changes = defaultdict(lambda: defaultdict(int))
        # and by the base change for adducts
        adduct_base_changes = defaultdict(lambda: defaultdict(int))

        # Iterate over each SSC
        for family_id, dsc in self.dsc_info.items():

            # If the family is not in the list to keep
            if family_id not in keep_families:

                # Skip it
                continue

            # Increment the number of reads
            total_counts['ssc'] += 1
            chr_counts[dsc['ref_name']]['ssc'] += 1

            # Increment the number of bases
            total_counts['bases'] += dsc['nbases']
            chr_counts[dsc['ref_name']]['bases'] += dsc['nbases']

            # Increment the number of adducts
            total_counts['adducts'] += len(dsc['adducts'])
            chr_counts[dsc['ref_name']]['adducts'] += len(dsc['adducts'])

            # Iterate over each of the positions with a SNP
            for snp_info in dsc['variants'].values():

                # Increment the number of SNPs
                total_counts['snps'] += 1
                chr_counts[dsc['ref_name']]['snps'] += 1

                # Increment the individual base change
                snp_base_changes[snp_info['var']][snp_info['ref']] += 1

            # Iterate over each of the positions with an adduct
            for adduct_info in dsc['adducts'].values():

                # Increment the number of adducts
                total_counts['adducts'] += 1
                chr_counts[dsc['ref_name']]['adducts'] += 1

                # Increment the individual base change
                adduct_base_changes[adduct_info['var']][adduct_info['ref']] += 1

        # Add all of the subset data to the totals
        total_counts['specimen'] = specimen
        total_counts['by_chr'] = chr_counts
        total_counts['snp_base_changes'] = snp_base_changes
        total_counts['adduct_base_changes'] = adduct_base_changes

        # Format the base change data as DataFrames
        snp_base_changes = pd.DataFrame(snp_base_changes).reindex(
            columns=['A', 'T', 'C', 'G'],
            index=['A', 'T', 'C', 'G']
        ).fillna(0).applymap(int)

        adduct_base_changes = pd.DataFrame(adduct_base_changes).reindex(
            columns=['A', 'T', 'C', 'G'],
            index=['A', 'T', 'C', 'G']
        ).fillna(0).applymap(int)

        return total_counts, pd.DataFrame(chr_counts), snp_base_changes, adduct_base_changes

    def write_total_json(self, folder=None, keep_families=None):
        """Save the total information to JSON."""

        fpo = os.path.join(folder, "total.json.gz")

        logger.info(f"Writing all output to {fpo}")

        with gzip.open(fpo, "wt") as handle:
            json.dump(
                {
                    family_id: dsc
                    for family_id, dsc in self.dsc_info.items()
                    if family_id in keep_families
                },
                handle
            )

    def write_adduct_family_list(self, folder=None, keep_families=None):
        """Save the list of all families which contain adducts."""

        # Output path for a text file containing the names of all families
        # which contain adducts
        fpo = os.path.join(folder, f"{specimen}.adduct.families.txt.gz")

        logger.info(f"Writing out {len(keep_families):,} families to {fpo}")
        with gzip.open(fpo, "wt") as handle:
            handle.write("\n".join(list(keep_families)))

    def write_adducts_gtf(self, folder=None, keep_families=None):
        """Write out the adduct information in GTF format."""

        # Format the output filepath
        fpo = os.path.join(folder, f"{specimen}.adduct.gtf")

        logger.info(f"Writing out adducts to {fpo} for {len(keep_families):,} families")

        # Make a list of all adducts
        adducts = [
            dict(
                seqname=dsc['ref_name'],
                # 0-index -> 1-index
                start=adduct_pos + 1,
                end=adduct_pos + 1,
                strand="+" if adduct_info["strand"] == "pos" else "-",
                attribute='adduct "%s"; read_as "%s"' % (adduct_info['ref'], adduct_info['var'])
            )
            for family_id, dsc in self.dsc_info.items()
            if family_id in keep_families
            for adduct_pos, adduct_info in dsc["adducts"].items()
        ]

        # If there are no adducts
        if len(adducts) == 0:

            print("No adducts found, skipping")
            return

        else:

            print(f"Writing out {len(adducts):,} adducts in GTF format")

        # Convert to a DataFrame
        adducts = pd.DataFrame(adducts)

        # Add in the fixed columns and resort
        adducts = adducts.assign(
            source=specimen,
            feature="adduct",
            score='.',
            frame='.'
        ).reindex(
            columns=[
                "seqname",
                "source",
                "feature",
                "start",
                "end",
                "score",
                "strand",
                "frame",
                "attribute"
            ]
        ).sort_values(
            by=['seqname', 'start']
        )

        # Write out as TSV
        adducts.to_csv(
            fpo,
            sep="\t",
            index=None,
            quoting=3
        )

ParseSSC(specimen)