#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Perform quality trimming on the input FASTQ data
process fixed_trim {
    container "${params.container__cutadapt}"
    cpus 1
    
    input:
    tuple val(specimen), path(R1), path(R2)

    output:
    tuple val(specimen), path("${R1.name.replaceAll(/.fastq.gz/, '')}.trimmed.fastq.gz"), path("${R2.name.replaceAll(/.fastq.gz/, '')}.trimmed.fastq.gz"), emit: reads
    tuple val(specimen), path("${specimen}.cutadapt.json"), emit: log

    script:
    template 'fixed_trim.sh'

}

// Combine all FASTQC data into a single report
process multiqc {
    container "${params.container__multiqc}"
    publishDir "${params.output}/3_end_trimmed/fastqc/", mode: 'copy', overwrite: true
    
    input:
    path "*"

    output:
    path "multiqc_report.html"

    script:
    template 'multiqc.sh'

}

// Assess quality of reads after trimming a fixed length
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

workflow trim_wf{

    take:
    reads_ch
    // tuple val(specimen), path(read_1), path(read_2)

    main:

    // Trim a fixed number of bases from the 5' end
    fixed_trim(reads_ch)

    // Run FASTQC on the reads post-trimming
    fastqc(
        fixed_trim.out.reads
    )

    // Assemble a multiqc report on the FASTQC data
    multiqc(
        fastqc.out.zip.flatten().toSortedList()
    )

    emit:
    reads = fixed_trim.out.reads

}