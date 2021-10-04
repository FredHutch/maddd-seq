#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Set default parameters
params.help = false
params.manifest = false

params.barcode_length = 5
params.minimum_alignment_score = 20
params.repeat_masker = false
params.bed = false
params.max_barcode_mismatch = 2
params.max_family_offset = 5
params.min_reads_per_ssc = 3
params.trim_length = 5

// Set the containers to use for each component
params.container__cutadapt = "quay.io/biocontainers/cutadapt:3.5--py36hc5360cc_0"
params.container__fastqc = "quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1"

// Import sub-workflows
include { manifest_wf } from './modules/manifest'
include { quality_wf } from './modules/quality'
include { barcodes_wf } from './modules/barcodes'
include { trim_wf } from './modules/trim'
include { align_wf } from './modules/align'
include { family_wf } from './modules/family'
include { variants_wf } from './modules/variants'

// Main workflow
workflow {
    
    // Parse the manifest
    manifest_wf(
        params.sample_sheet
    )
    // output:
    //   reads:
    //     tuple val(specimen), path(read_1), path(read_2)

    // Perform quality trimming on the ends of all reads
    quality_wf(
        manifest.out.reads
    )
    // output:
    //   reads:
    //     tuple val(specimen), path(read_1), path(read_2)
    // publish:
    //   1_input_data/input.fastqc.html
    //   1_input_data/quality_trimmed.fastqc.html

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
    //   2_barcode_trimmed/barcode_frequency.csv

    // Trim a fixed number of bases from the beginning of each read
    trim_wf(
        barcodes_wf.out.reads
    )
    // output:
    //   bam:
    //     tuple val(specimen), path(bam)
    // publish:
    //   3_end_trimmed/end_trimmed.fastqc.html

    // Align the barcode-clipped reads to the reference genome
    align_wf(
        trim_wf.out.bam
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

}