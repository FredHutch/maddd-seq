#!/usr/bin/env python3
"""
Parse the SSC data in order to:
  - Construct DSCs BAMs which combine both strands
  - Call adducts from the SSC data
  - Call SNPs from the SSC data
  - Summarize the total number of adducts and SNPs
  - Summarize the number of adducts and SNPs per chromosome
  - Summarize the number of adducts and SNPs per position within each read
"""

from collections import defaultdict
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


# Output path for the total and summarized data
total_output = "total.json.gz"
print(f"Output path: {total_output}")
summary_output = "summary.json.gz"
print(f"Output path: {summary_output}")

# Output path for the data grouped by contig
chr_output = "by_chr.csv.gz"
print(f"Output path: {chr_output}")

# Output path for SNP data grouped by base change
snps_base_output = "snps_by_base.csv.gz"
print(f"Output path: {snps_base_output}")

# Output path for adduct data grouped by base change
adduct_base_output = "adducts_by_base.csv.gz"
print(f"Output path: {adduct_base_output}")

# Output path for variant data organized by position in the reads
base_positions_output = "by_read_position.csv.gz"
print(f"Output path: {base_positions_output}")

def parse_read(read):
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
        # Variants and reference positions
        **parse_variants(read)
    )

    # Tuple, Family ID, orientation, and a dict with position and variants
    return family_id, orient, read_details


def parse_variants(read):
    """
    For any position in which the read differs from the reference,
    include the reference position and the variant base in a dict.
    """

    # Encode the output as two dicts, one for variants and one for references
    variants = dict()
    references = dict()

    # Get the reference sequence
    rseq = read.get_reference_sequence()

    # get_aligned_pairs() returns a tuple of positions
    for qpos, rpos in read.get_aligned_pairs():

        # If there is an indel
        if rpos is None or qpos is None:
            
            # Skip it
            continue

        # Get the aligned and reference bases
        qbase = read.query_sequence[qpos]
        rbase = rseq[rpos - read.reference_start]

        # If the query matches the reference
        if qbase == rbase:

            # Skip it
            continue

        # Otherwise, record the query base in the output
        variants[rpos] = qbase
        # As well as the reference position at that location
        references[rpos] = rbase

    return dict(
        variants=variants,
        references=references
    )


def parse_bam(bam_fp):
    """Parse all of the data from a single BAM."""

    # Populate the outputs as a dict of dicts
    # First key is the family ID, second key is the orientation
    ssc_data = defaultdict(dict)

    # Open the input BAM file for reading
    with pysam.AlignmentFile(bam_fp, "rb") as bam:

        # For each of the reads on the positive strand    
        for read in bam:

            # Get the stats for this read
            family_id, orientation, read_stats = parse_read(read)

            # Add to the dictionary
            ssc_data[family_id][orientation] = read_stats

    return ssc_data

def merge_strands(pos_ssc, neg_ssc):
    """Merge information from both strands to identify SNPs and adducts."""

    # Output data is a dictionary keyed by family
    output = dict()

    # Also keep track of the number of mutations and adducts
    # as a function of the position in the read
    base_positions = dict(
        nreads=defaultdict(int),
        snps=defaultdict(int),
        adducts=defaultdict(int)
    )

    # Iterate over each family
    for family_id in pos_ssc:

        # Make sure that the family is in the negative strand data as well
        assert neg_ssc.get(family_id) is not None, f"Missing negative strand information for {family_id}"

        # Merge data from fwd/rev, pos/neg
        output[family_id] = merge_single_family(pos_ssc.get(family_id), neg_ssc.get(family_id))

        # Add the nreads/snps/adducts data as a function of read position
        add_base_position_info(
            output[family_id],
            pos_ssc.get(family_id),
            neg_ssc.get(family_id),
            base_positions
        )

        # Add the number of aligned bases
        output[family_id]['nbases'] = calc_nbases(pos_ssc.get(family_id), neg_ssc.get(family_id))

    return output, base_positions


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

        # If both strands agree
        if pos_base == neg_base:

            # Then this is a mutation
            output['mutations'][var_pos] = pos_base
            output['references'][var_pos] = reference_bases[var_pos]

        # if the strands do not agree, and the positive strand is unchanged
        elif pos_base == 'ref':
            # Add the adduct to the list
            output['adducts'][var_pos] = neg_base, 'neg'
            output['references'][var_pos] = reference_bases[var_pos]

        # If the negative strand is unchanged
        elif neg_base == 'ref':
            # Add the adduct to the list
            output['adducts'][var_pos] = pos_base, 'pos'
            output['references'][var_pos] = reference_bases[var_pos]

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


def format_output(ssc_dat):
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
    for _, family_dat in ssc_dat.items():

        # Increment the number of reads
        total_counts['ssc'] += 1
        chr_counts[family_dat['ref_name']]['ssc'] += 1

        # Increment the number of bases
        total_counts['bases'] += family_dat['nbases']
        chr_counts[family_dat['ref_name']]['bases'] += family_dat['nbases']

        # Increment the number of adducts
        total_counts['adducts'] += len(family_dat['adducts'])
        chr_counts[family_dat['ref_name']]['adducts'] += len(family_dat['adducts'])

        # Iterate over each of the positions with a SNP
        for snp_pos, snp_base in family_dat['mutations'].items():

            # Increment the number of SNPs
            total_counts['snps'] += 1
            chr_counts[family_dat['ref_name']]['snps'] += 1

            # Increment the individual base change
            snp_base_changes[family_dat['references'][snp_pos]][snp_base] += 1

        # Iterate over each of the positions with an adduct
        for adduct_pos, (adduct_base, adduct_strand) in family_dat['adducts'].items():

            # Increment the number of adducts
            total_counts['adducts'] += 1
            chr_counts[family_dat['ref_name']]['adducts'] += 1

            # Get the base of the reference at this position
            ref_base = family_dat['references'][adduct_pos]

            # If the adduct is on the positive strand
            if adduct_strand == 'pos':

                # Increment the individual base change
                adduct_base_changes[ref_base][adduct_base] += 1

            # If the adduct is on the negative strand
            else:

                # Increment the individual base change from the other strand
                adduct_base_changes[
                    dict(
                        A='T',
                        T='A',
                        C='G',
                        G='C'
                    )[ref_base]
                ][adduct_base] += 1

    # Add all of the subset data to the totals
    total_counts['specimen'] = "${specimen}"
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

# Parse the information from the positive strand
pos_ssc = parse_bam(input_pos_bam)
# Parse the information from the negative strand
neg_ssc = parse_bam(input_neg_bam)

# Join all of the information from both strands
ssc_dat, base_positions = merge_strands(pos_ssc, neg_ssc)

# Save the total information to JSON
with gzip.open(total_output, "wt") as handle:
    json.dump(ssc_dat, handle)

# Format all of the output, both by contig and overall
summary_dat, by_chr, snp_base_changes, adduct_base_changes = format_output(ssc_dat)

# Save the summary information to JSON
with gzip.open(summary_output, "wt") as handle:
    json.dump(summary_dat, handle)

# Save the by-chr information to CSV
by_chr.to_csv(chr_output)

# Save the SNP by-base information to CSV
snp_base_changes.to_csv(snps_base_output)

# Save the adduct by-base information to CSV
adduct_base_changes.to_csv(adduct_base_output)

# Format the data by read position as a DataFrame
base_positions = pd.DataFrame(base_positions).fillna(0).applymap(int)
base_positions.to_csv(base_positions_output)
