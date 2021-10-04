#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Assess quality of input data
process fastqc_input {
    container "${params.container__fastqc}"
    cpus 1
    
    input:
    tuple val(specimen), path(R1), path(R2)

    output:
    tuple val(specimen), path("fastqc/fastqc.html")

    script:
    template 'fastqc.sh'

}

// Assess quality of trimmed data
process fastqc_trimmed {
    container "${params.container__fastqc}"
    cpus 1
    
    input:
    tuple val(specimen), path(R1), path(R2)

    output:
    tuple val(specimen), path("fastqc/fastqc.html")

    script:
    template 'fastqc.sh'

}

// Perform quality trimming on the input FASTQ data
process quality_trim {
    container "${params.container__cutadapt}"
    cpus 1
    
    input:
    tuple val(specimen), path(R1), path(R2)

    output:
    tuple val(specimen), path("trimmed.R1.fastq.gz"), path("trimmed.R1.fastq.gz"), emit: reads
    tuple val(specimen), path("${specimen}.cutadapt.json"), emit: log

    script:
    template 'quality_trim.sh'

}

workflow quality_wf{

    take:
    reads_ch
    // tuple val(specimen), path(read_1), path(read_2)

    main:

    // Generate quality metrics for the input data
    fastqc_input(reads_ch)

    // Run quality trimming
    quality_trim(reads_ch)

    // Generate quality metrics for the trimmed data
    fastqc_trim(reads_ch)

    emit:
    reads = quality_trim.out.reads

}