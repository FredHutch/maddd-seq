# Mutation and DNA Damage Detection by massively-parallel sequencing (MADDD-seq)

## Background

When analyzing the variants present in a collection of nucleic acids, one option is
to perform "MADDD-seq" in which both strands of the DNA duplex are sequenced
independently. This method for preparing nucleic acids for high-throughput sequencing
relies on an approach in which unique barcodes are attached to both ends of each
distinct piece of DNA. By grouping together sequence reads which share the same
barcode combinations, the analysis can generate highly accurate measurements of the
nucleotide sequence of each independent strand. Using that strand-specific information,
the presence of single-nucleotide polymorphisms (SNPs) and adducts is inferred for
all of the DNA fragments which were analyzed.

One of the primary determinents of the quality of this type of sequencing
analysis is the combination of filtering parameters which is applied to the raw sequence
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
3. Trim an additional N bases from the beginning of each read to reduce errors introduced by overhang filling
4. Align the barcode-trimmed reads to the reference genome
5. Group together sequences which share the same barcode and which are aligned to the same position into 'families'
6. Infer the Single-Strand Consensus (SSC) from each family of input data
7. Compare the SSC sequences from each strand to generate the double-strand consensus (DSC) sequence
8. Summarize the number and position of DNA variants and adducts detected

The options available for this analysis are:

- Quality threshold used for trimming in (1)
- Length of the barcodes trimmed in (2)
- List of barcode sequences used for error-correction in (2)
- The maximum number of mismatches used for barcode error-correction in (2)
- Additional number of bases to remove from the beginning of each read in (3)
- Minimum threshold for alignment scores in (4)
- A set of genomic regions used for target capture (used to limit the regions of analysis in 4)
- The maximum distance between the alignment start position of reads grouped together in (5)
- Minimum number of reads per SSC used in analysis after (6)
- An optional set of genomic coordinates to omit while detecting SNPs and adducts in (8)

## Outputs

```
1_input_data/input/multiqc_report.html
Summary of sequence quality for input data

1_input_data/quality_trimmed/multiqc_report.html
Summary of sequence quality after trimming by quality score

2_barcode_trimmed/{specimen}/barcode_counts.csv.gz
Summary of the number of reads from each specimen with each uncorrected barcode

2_barcode_trimmed/{specimen}/barcode_corrections.csv.gz
Summary of the number of reads from each specimen with each corrected barcode
(after barcode error correction)

2_barcode_trimmed/{specimen}/{specimen}.barcodes.pdf
Figure displaying the distribution of the number of reads sequenced per barcode

2_barcode_trimmed/multiqc_report.html
Summary of sequence quality after trimming barcodes

3_end_trimmed/multiqc_report.html
Summary of sequence quality after trimming additional bases from the 5' end of each read

4_aligned/multiqc_report.html
Summary of the overall number of reads aligned per specimen

4_aligned/{specimen}/{specimen}.flagstats
Summary of the number of reads aligned per contig per specimen

6_all_SSC/<specimen>/[NEG|POS].SSC.bam[.bai]
Aligned SSCs (unfiltered by number of reads) per specimen, both positive and negative strands

6_all_SSC/<specimen>/SSC.details.csv.gz
For each SSC (unfiltered by number of reads) in a specimen, details including the number of reads

7_filtered_SSC/<specimen>/alignments/DSC.bam[.bai]
Aligned DSCs per specimen, after filtering by the minimum number of reads per SSC

7_filtered_SSC/<specimen>/alignments/[POS|NEG].SSC.bam[.bai]
Aligned SSCs per specimen, after filtering by the minimum number of reads per SSC

7_filtered_SSC/<specimen>/alignments/DSC.vcf.gz
Tabular summary of variant positions detected per specimen

7_filtered_SSC/<specimen>/stats/{specimen}.snps_by_base.csv.gz
Number of SNPs detected as a function of what bases were changed (e.g. A -> C)

7_filtered_SSC/<specimen>/stats/{specimen}.adducts_by_base.csv.gz
Number of adducts detected as a function of what bases were changed (e.g. A -> C)

7_filtered_SSC/<specimen>/stats/{specimen}.by_chr.csv.gz
Number of SNPs and adducts detected per chromosome (contig)

7_filtered_SSC/<specimen>/stats/{specimen}.by_read_position.csv.gz
Number of SNPs and adducts detected as a function of the position along the read

7_filtered_SSC/<specimen>/stats/{specimen}.SSC.csv.gz
Tabular summary of all filtered SSCs, including SNPs, adducts, and read length

7_filtered_SSC/plots/DSC.summary.pdf
Visual summary of various metrics summarizing the filtered DSCs

7_filtered_SSC/plots/DSC.summary.pdf
Visual summary of the distribution of SNPs and adducts as a function of position along the read
```
