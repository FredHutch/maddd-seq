#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Remove the UMIs from each read pair and move to the header
process clip_barcodes {
    container "${params.container__pysam}"

    input:
    tuple val(specimen), path(R1), path(R2)

    output:
    tuple val(specimen), path("${R1.name.replaceAll(/.fastq.gz/, '')}.clipped.fastq.gz"), path("${R2.name.replaceAll(/.fastq.gz/, '')}.clipped.fastq.gz")

    script:
    template 'clip_barcodes.py'

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

workflow barcodes_wf{

    take:
    fastq_ch

    main:

    // Remove the barcodes from each end of the read
    // and move to the header of the FASTQ
    clip_barcodes(fastq_ch)

    // Assess the quality of the WGS data after barcodes are removed
    fastqc(
        clip_barcodes.out
    )

    // Combine all of the FASTQC reports into a single file
    multiqc(
        fastqc.out.zip
    )

    // // Count up the number of reads with each barcode
    // count_barcodes(clip_barcodes.out)

    // // Perform error correction on the barcodes
    // correct_barcode_errors(count_barcodes.out)

    // // Modify the barcode in the unaligned BAM to
    // // the corrected sequence
    // update_barcodes(
    //     fastq_ch.combine(
    //         correct_barcode_errors.out
    //     )
    // )

    emit:
    bam = fastq_ch

}