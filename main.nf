#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Set default parameters
params.help = false
params.sample_sheet = false
params.output = false
params.genome = false

// Quality trimming
params.min_qvalue = 20
params.min_align_score = 40

// Unique molecular tags
params.barcode_length = 12
params.barcode_max_homopolymer = 6
params.max_barcode_mismatch = 2

// Trim a fixed amount from the 5' of both reads
params.trim_length = 5

params.repeat_masker = false
params.bed = false
params.max_family_offset = 5
params.min_reads_per_ssc = 3

// Set the containers to use for each component
params.container__cutadapt = "quay.io/biocontainers/cutadapt:3.5--py36hc5360cc_0"
params.container__fastqc = "quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1"
params.container__multiqc = "quay.io/biocontainers/multiqc:1.11--pyhdfd78af_0"
params.container__pysam = "quay.io/biocontainers/pysam:0.17.0--py36h61e5637_0"
params.container__pandas = "quay.io/fhcrc-microbiome/python-pandas:4a6179f"
params.container__python_plotting = "quay.io/hdc-workflows/python-plotting:b50a842"
params.container__bwa = "quay.io/hdc-workflows/bwa-samtools:latest"

// Import sub-workflows
include { manifest_wf } from './modules/manifest'
include { quality_wf } from './modules/quality'
include { barcodes_wf } from './modules/barcodes'
include { trim_wf } from './modules/trim'
include { align_wf } from './modules/align'
include { family_wf } from './modules/family'
include { variants_wf } from './modules/variants'


// Function which prints help message text
def helpMessage() {
    log.info"""
Usage:

nextflow run FredHutch/wgs-duplex-seq <ARGUMENTS>

Required Arguments:
  --sample_sheet        CSV file listing samples with headers: specimen, R1, and R2
  --output              Folder to write output files to
  --genome              Reference genome indexed for alignment by BWA

Optional Arguments:
  --min_qvalue          Minimum quality score used to trim data (default: ${params.min_qvalue})
  --min_align_score     Minimum alignment score (default: ${params.min_align_score})
  --barcode_length      Length of the barcodes ligated to the 5' end of each read
                        (default: ${params.barcode_length})
  --barcode_max_homopolymer
                        Maximum number of repeated bases in barcode sequences
                        (default: ${params.barcode_max_homopolymer})
  --max_barcode_mismatch
                        Maximum number of nucleotide differences used for barcode
                        error correction
  --trim_length         Number of bases to trim from the 5' of each read after
                        ligated barcode sequences are removed

Manifest:
  The manifest is a CSV listing all of the duplex sequencing data to be analyzed.
  The manifest must contain the column headers: specimen,R1,R2
  The R1 and R2 file paths may start with ftp://, s3://, or even just be a path to a local file.
  The specimen must be unique, and only contain a-z, A-Z, 0-9, or _.
    """.stripIndent()
}


// Main workflow
workflow {

    // Show help message if the user specifies the --help flag at runtime
    // or if --sample_sheet and --output are not provided
    if ( params.help || params.sample_sheet == false || params.output == false || params.genome == false ){
        // Invoke the function above which prints the help message
        helpMessage()
        // Exit out and do not run anything else
        exit 1
    }

    // Parse the manifest
    manifest_wf(
        params.sample_sheet
    )
    // output:
    //   reads:
    //     tuple val(specimen), path(read_1), path(read_2)

    // Perform quality trimming on the ends of all reads
    quality_wf(
        manifest_wf.out.reads
    )
    // output:
    //   reads:
    //     tuple val(specimen), path(read_1), path(read_2)
    // publish:
    //   1_input_data/input/multiqc_report.html
    //   1_input_data/quality_trimmed/multiqc_report.html

    // Clip the barcodes from the ends of the FASTQ files
    // This part of the workflow will also perform error
    // correction on the barcode sequences.
    barcodes_wf(
        quality_wf.out.reads
    )
    // output:
    //   bam:
    //     tuple val(specimen), path(bam)
    // publish:
    //   2_barcode_trimmed/{specimen}/barcode_counts.csv.gz
    //   2_barcode_trimmed/{specimen}/barcode_corrections.csv.gz
    //   2_barcode_trimmed/multiqc_report.html

    // Trim a fixed number of bases from the beginning of each read
    trim_wf(
        barcodes_wf.out.reads
    )
    // output:
    //   bam:
    //     tuple val(specimen), path(bam)
    // publish:
    //   3_end_trimmed/fastqc/multiqc_report.html
    //   3_end_trimmed/cutadapt/multiqc_report.html

    // Align the barcode-clipped reads to the reference genome
    align_wf(
        trim_wf.out.reads,
        Channel.fromPath("${params.genome}").toSortedList()
    )
    // output:
    //   bam:
    //     tuple val(specimen), path(bam)
    // publish:
    //   4_aligned/alignment_summary.csv

    // Group reads into families based on barcodes and alignment position
    // This sub-workflow will also collapse and summarize SSCs and DSCs
    family_wf(
        align_wf.out.bam
    )
    // output:
    //   bam: 
    //     tuple val(specimen), path(ssc_bam), path(dsc_bam)
    // publish:
    //   5_families/family_summary.csv
    //   5_families/<specimen>/aligned.bam[.bai]
    //   6_SSC/<specimen>/SSC.bam[.bai]
    //   6_SSC/<specimen>/SSC.details.csv.gz
    //   7_DSC/<specimen>/DSC.bam[.bai]
    //   7_DSC/<specimen>/DSC.details.csv.gz
    //   7_DSC/<specimen>/DSC.positional.csv.gz

    // Call variants and adducts
    variants_wf(
        family_wf.out.bam
    )
    // publish:
    //   8_variants/<specimen>/variants.vcf
    //   8_variants/<specimen>/adducts.vcf
    //   8_variants/variant_summary.csv
    //   2_barcode_trimmed/multiqc_report.html
    //   3_end_trimmed/fastqc/multiqc_report.html
    //   3_end_trimmed/cutadapt/multiqc_report.html

}