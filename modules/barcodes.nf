#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Remove the UMIs from each read pair and move to the header
process clip_barcodes {
    container "${params.container__cutadapt}"
    publishDir "${params.output}/2_barcode_trimmed/", mode: 'copy', overwrite: true, pattern: "barcode_counts.csv.gz"

    input:
    tuple val(specimen), path(R1), path(R2)

    output:
    tuple val(specimen), path("${R1.name.replaceAll(/.fastq.gz/, '')}.clipped.fastq.gz"), path("${R2.name.replaceAll(/.fastq.gz/, '')}.clipped.fastq.gz"), emit: reads
    tuple val(specimen), path("barcode_counts.csv.gz"), emit: counts
    tuple val(specimen), path("${specimen}.cutadapt.json"), emit: json

    script:
    template 'clip_barcodes.sh'

}

// Assess quality of reads after removing barcodes
process fastqc {
    container "${params.container__fastqc}"
    
    input:
    tuple val(specimen), path(R1), path(R2)

    output:
    path "fastqc/*.zip", emit: zip
    path "fastqc/*.html", emit: html

    script:
    template 'fastqc.sh'

}

// Combine all FASTQC data into a single report
process multiqc {
    container "${params.container__multiqc}"
    publishDir "${params.output}/2_barcode_trimmed/", mode: 'copy', overwrite: true
    
    input:
    path "*"

    output:
    path "multiqc_report.html"

    script:
    template 'multiqc.sh'

}

// Perform error correction on the barcodes
process correct_barcode_errors {
    container "${params.container__pandas}"
    publishDir "${params.output}/2_barcode_trimmed/", mode: 'copy', overwrite: true, glob: "barcode_corrections.csv.gz"

    input:
    tuple val(specimen), path("barcode_counts.csv.gz")

    output:
    tuple val(specimen), path("barcode_corrections.csv.gz")

    script:
    template 'correct_barcode_errors.py'

}

// Apply the corrected barcode sequences to the FASTQ data
process update_barcodes {
    container "${params.container__pysam}"

    input:
    tuple val(specimen), path(R1), path(R2), path("barcode_corrections.csv.gz")

    output:
    tuple val(specimen), path("${R1.name.replaceAll(/.fastq.gz/, '')}.corrected.fastq.gz"), path("${R2.name.replaceAll(/.fastq.gz/, '')}.corrected.fastq.gz"), emit: reads

    script:
    template 'update_barcodes.py'

}

workflow barcodes_wf{

    take:
    fastq_ch

    main:

    // Remove the barcodes from each end of the read
    // and move to the header of the FASTQ
    clip_barcodes(fastq_ch)

    // Assess the quality of the WGS data after barcodes are removed
    fastqc(
        clip_barcodes.out.reads
    )

    // Combine all of the FASTQC reports into a single file
    multiqc(
        fastqc.out.zip.flatten().toSortedList()
    )

    // Perform error correction on the barcodes
    correct_barcode_errors(
        clip_barcodes.out.counts
    )

    // // Modify the barcode in the FASTQ files to
    // // the corrected sequence
    // update_barcodes(
    //     fastq_ch.combine(
    //         correct_barcode_errors.out
    //     )
    // )

    emit:
    bam = fastq_ch

}