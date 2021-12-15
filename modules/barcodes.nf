#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Remove the UMIs from each read pair and move to the header
process clip_barcodes {
    container "${params.container__cutadapt}"
    publishDir "${params.output}/2_barcode_trimmed/${specimen}/", mode: 'copy', overwrite: true, pattern: "barcode_counts.csv.gz"
    publishDir "${params.output}/2_barcode_trimmed/${specimen}/clip_barcodes_intermediate/", mode: 'copy', pattern: "*", enabled: params.save_intermediates
    label "cpu_medium"

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
    publishDir "${params.output}/2_barcode_trimmed/${specimen}/fastqc_intermediate/", mode: 'copy', enabled: params.save_intermediates
    label "io_limited"
    
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
    label "io_limited"
    
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
    publishDir "${params.output}/2_barcode_trimmed/${specimen}/", mode: 'copy', overwrite: true, glob: "barcode_corrections.csv.gz"
    label "mem_medium"

    input:
    tuple val(specimen), path("barcode_counts.csv.gz")
    path "barcodes.txt"

    output:
    tuple val(specimen), path("barcode_corrections.csv.gz")

    script:
    template 'correct_barcode_errors.py'

}

// Make plots for the barcodes, both raw and corrected
process plot_barcodes {
    container "${params.container__python_plotting}"
    publishDir "${params.output}/2_barcode_trimmed/${specimen}/", mode: 'copy', overwrite: true
    label "io_limited"

    input:
    tuple val(specimen), path("barcode_counts.csv.gz"), path("barcode_corrections.csv.gz")

    output:
    path "${specimen}.barcodes.pdf", emit:pdf
    path "${specimen}.corrected_barcodes.csv", emit:csv

    script:
    template 'plot_barcodes.py'

}


workflow barcodes_wf{

    take:
    fastq_ch
    barcodes_txt

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
        clip_barcodes.out.counts,
        barcodes_txt
    )

    // Make plots summarizing the barcodes
    plot_barcodes(
        clip_barcodes.out.counts.join(
            correct_barcode_errors.out
        )
    )

    emit:
    // The reads have not yet had their barcodes corrected, but
    // they have been clipped from the sequence and moved to a header
    reads = clip_barcodes.out.reads
    // The CSV maps each barcode to its corrected form
    csv = correct_barcode_errors.out
    // The channel contains 1 CSV per specimen, which contains the number of reads per corrected barcode
    counts = plot_barcodes.out.csv

}