# WGS Duplex Seq
Analysis of whole-genome shotgun duplex sequencing data

## Background

When analyzing the variants present in a collection of nucleic acids, one option is
to perform "duplex sequencing" in which both strands of the DNA duplex are sequenced
independently. This method for preparing nucleic acids for high-throughput sequencing
relies on an approach in which unique barcodes are attached to both ends of each
distinct piece of DNA. By grouping together sequence reads which share the same
barcode combinations, the analysis can generate highly accurate measurements of the
nucleotide sequence of each independent strand.

One of the primary determinents of the quality of this type of duplex sequencing
analysis is the combination of filtering parameters which is applied to the raw
data. The elements of filtering which are most important are (a) the minimum number
of reads required per DNA strand and (b) the number of bases which are trimmed from
each end of the DNA molecule prior to variant calling. The underlying reason for the
end-trimming is that the barcode ligation process involves a gap repair process which
is highly error prone. By removing the sequences which are generated from the end of
the DNA molecule, those gap-filled bases can be excluded from downstream analysis
as is appropriate.

After identifying the consensus sequence from each strand of each molecule, the
strands can be compared against each other to identify positions at which there is
a mismatch within the duplex molecule. Those double-stranded mismatches reflect
DNA adducts which were originally present in the genomic material which was processed.
One of the goals of this analysis pipeline is to identify those adducts, both in
terms of the total rate of adducts identified in a particular specimen, as well as
the position of those adducts in the genome used for alignment.

## Approach

The steps of this analysis workflow are:

1. Perform quality trimming and QC of the input whole-genome shotgun sequence (WGS) data
2. Trim barcodes from the ends of WGS reads and tag each pair with the concatenated barcode
3. Trim an additional N bases from the beginning of each read to reduce errors introduced by blunt-end ligation
4. Align the barcode-trimmed reads to the reference genome
5. Group together sequences which share the same barcode and which are aligned to the same position into 'families'
6. Infer the Single-Strand Consensus (SSC) from each family of input data
7. Compare the SSC sequences from each strand to generate the double-strand consensus (DSC) sequence
8. Trim the ends from each DSC and summarize the number and position of DNA variants and adducts detected

The options available for this analysis are:

- Quality threshold used for trimming in (1)
- Length of the barcodes trimmed in (2)
- Additional number of bases to remove from the beginning of each read in (3)
- Minimum threshold for alignment scores in (4)
- Whether or not repetitive elements should be masked in (4)
- A set of genomic regions used for target capture (used to limit the regions of analysis in 4)
- The maximum number of mismatches used for grouping molecular barcodes in (5)
- The maximum distance between the alignment start position of reads grouped together in (5)
- Minimum number of reads per SSC used in analysis after (6)
- Minimum number of DSCs required to call variants in (8)

## Outputs

```
1_input_data/input.fastqc.html
Summary of sequence quality for input data

1_input_data/quality_trimmed.fastqc.html
Summary of sequence quality after trimming by quality score

2_barcode_trimmed/barcode_frequency.csv
Summary of the number of reads from each specimen with each barcode
(after barcode error correction)

3_end_trimmed/end_trimmed.fastqc.html
Summary of sequence quality after trimming additional 5' bases

4_aligned/alignment_summary.csv
Summary of the number of reads aligned per contig per specimen

5_families/family_summary.csv
Summary of the distribution of family sizes per specimen

5_families/<specimen>/aligned.bam[.bai]
Aligned reads, tagged by barcode and family

6_SSC/<specimen>/SSC.bam[.bai]
Aligned SSCs per specimen

6_SSC/<specimen>/SSC.details.csv.gz
For each SSC in a specimen, details including the number of reads and mismatches

7_DSC/<specimen>/DSC.bam[.bai]
Aligned DSCs per specimen

7_DSC/<specimen>/DSC.details.csv.gz
For each DSC in a specimen, details including the number of reads on each strand, mismatches, and adducts

7_DSC/<specimen>/DSC.positional.csv.gz
Over all DSCs in a specimen, the rate of variants and adducts as a function of read position

8_variants/<specimen>/variants.vcf
Position of variants detected per specimen

8_variants/<specimen>/adducts.vcf
Position of adducts detected per specimen

8_variants/variant_summary.csv
Number and rate of variants and adducts per specimen
```
