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
4. Align the barcode-trimmed reads to the reference genome and trim any 3' bases which overhang the paired end
5. Infer the Single-Strand Consensus (SSC) from each family of input data which share the same barcode and which are aligned to the same position
6. Filter SSC data by sequencing depth and combine to generate the double-strand consensus (DSC) sequence, yielding the DNA variants and adducts detected

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
- An optional set of genomic coordinates to omit while detecting SNPs and adducts in (6)

## Outputs

```
├── 1_input_data
│   ├── input
│   │   |   # Summary of sequence quality for input data
│   │   └── multiqc_report.html
│   └── quality_trimmed
│       ├── <SPECIMEN>
│       │   |   # Quality trimmed reads
│       │   ├── <SPECIMEN>_R1_001.trimmed.fastq.gz
│       │   └── <SPECIMEN>_R2_001.trimmed.fastq.gz
│       │   # Summary of sequence quality after trimming by quality score
│       └── multiqc_report.html
├── 2_barcode_trimmed
│   ├── <SPECIMEN>
│   │   |   # Visual summary of the number of barcodes pre- and post-trimming
│   │   ├── <SPECIMEN>.barcodes.pdf
│   │   |   # Summary of the number of reads from each specimen with each corrected barcode
│   │   ├── barcode_corrections.csv.gz
│   │   |   # Summary of the number of reads from each specimen with each uncorrected barcode
│   │   └── barcode_counts.csv.gz
│   │   # Summary of sequence quality after trimming barcodes
│   └── multiqc_report.html
├── 3_end_trimmed
│   ├── <SPECIMEN>
│   │   |   # End-trimmed reads
│   │   ├── <SPECIMEN>_R1_001.trimmed.clipped.trimmed.fastq.gz
│   │   └── <SPECIMEN>_R2_001.trimmed.clipped.trimmed.fastq.gz
│   |   # Summary of sequence quality after 5' fixed-length end-trimming
│   └── multiqc_report.html
├── 4_aligned
│   |   # Summary of the overall number of reads aligned per specimen
│   ├── multiqc_report.html
│   ├── reads
│   │   |   # Aligned reads
│   │   ├── <SPECIMEN>.aligned.bam
│   │   ├── <SPECIMEN>.aligned.bam.bai
│   │   └── <SPECIMEN>.aligned.bam.csi
│   └── trim_overhang
│       └── <SPECIMEN>
│           |   # Overhang-trimmed reads
│           ├── <SPECIMEN>_R1.fastq.gz
│           └── <SPECIMEN>_R2.fastq.gz
├── 5_all_SSC
│   └── <SPECIMEN>
│       |   # Summary metrics for every family of reads
│       ├── <SPECIMEN>.unfiltered.SSC.details.csv.gz
│       |   # Aligned SSCs, negative strand
│       ├── NEG.SSC.bam
│       ├── NEG.SSC.bam.bai
│       |   # Aligned SSCs, positive strand
│       ├── POS.SSC.bam
│       └── POS.SSC.bam.bai
└── 6_filtered_SSC
    ├── <SPECIMEN>
    │   |   # Location of mutations and the depth of coverage at each position
    │   ├── DSC.tsv.gz
    │   ├── alignments
    │   │   |   # Aligned DSC data (BAM format, per read)
    │   │   ├── DSC.bam
    │   │   ├── DSC.bam.bai
    │   │   |   # Aligned DSC data (PILEUP format, per position)
    │   │   ├── DSC.pileup.gz
    │   │   |   # Tabular summary of variant positions detected per specimen
    │   │   ├── DSC.vcf.gz
    │   │   |   # Alignments of SSC data after filtering by minimum read number
    │   │   ├── NEG.SSC.bam
    │   │   ├── NEG.SSC.bam.bai
    │   │   ├── POS.SSC.bam
    │   │   └── POS.SSC.bam.bai
    │   └── stats
    │       |   # Location of mutations and the depth of coverage at each position
    │       ├── <SPECIMEN>.DSC.tsv.gz
    │       |   # Position, read length, and number of reads per filtered SSC
    │       ├── <SPECIMEN>.SSC.csv.gz
    │       |   # Adduct data with position and strand-specific mismatch information
    │       ├── <SPECIMEN>.adduct.gtf
    │       |   # Number of adducts detected as a function of what bases were changed (e.g. A -> C)
    │       ├── <SPECIMEN>.adducts_by_base.csv.gz
    │       |   # Number of SNPs and adducts detected per chromosome (contig)
    │       ├── <SPECIMEN>.by_chr.csv.gz
    │       |   # Number of SNPs and adducts detected as a function of the position along the read
    │       ├── <SPECIMEN>.by_read_position.csv.gz
    │       |   # Number of SNPs detected as a function of what bases were changed (e.g. A -> C)
    │       └── <SPECIMEN>.snps_by_base.csv.gz
    ├── plots
    │   |   # Collection of figures summarizing SNP and adduct rates across all specimens
    │   └── report.pdf
    ├── reads
    │   |   # Aligned reads from DSCs containing adducts
    │   ├── <SPECIMEN>.adduct.reads.bam
    │   ├── <SPECIMEN>.adduct.reads.bam.bai
    │   └── <SPECIMEN>.adduct.reads.bam.csi
    └── tables
        |   # Tabular summary of the number of adducts per specimen
        ├── adducts_by_base.csv
        |   # Tabular summary of the number of SNPs per specimen
        ├── snps_by_base.csv
        |   # Tabular summary of the SNP and adduct rate per specimen
        └── summary.csv
```
