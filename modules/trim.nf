#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Perform quality trimming on the input FASTQ data
process fixed_trim {
    container "${params.container__cutadapt}"
    cpus 1
    
    input:
    tuple val(specimen), path("trimmed.R1.fastq.gz"), path("trimmed.R1.fastq.gz"), emit: reads

    output:
    tuple val(specimen), path(R1), path(R2), emit: reads
    tuple val(specimen), path("${specimen}.cutadapt.json"), emit: log

    script:
    template 'fixed_trim.sh'

}

workflow trim_wf{

    take:
    reads_ch
    // tuple val(specimen), path(read_1), path(read_2)

    main:

    // Trim a fixed number of bases from the 5' end
    fixed_trim(reads_ch)

    emit:
    reads = fixed_trim.out.reads

}