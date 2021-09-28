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

1. Trim barcodes from the input whole-genome shotgun sequence (WGS) data
2. Align the barcode-trimmed reads to the reference genome
3. Group together sequences which share the same barcode and which are aligned to the same position into 'families'
4. Infer the Single-Strand Consensus (SSC) from each family of input data
5. Compare the SSC sequences from each strand to generate the double-strand consensus (DSC) sequence
6. Trim the ends from each DSC and summarize the number and position of DNA variants and adducts detected

The options available for this analysis are:

- Length of the barcodes trimmed in (1)
- Minimum threshold for alignment scores in (2)
- Whether repetitive elements should be masked in (2)
- A set of genomic regions used for target capture (used to limit the regions of analysis)
- The maximum number of mismatches used for grouping barcodes in (3)
- The maximum distance between the alignment start position of reads grouped together in (3)
- Minimum number of reads per SSC
- Number of bases to trim from the ends of each sequenced molecule
