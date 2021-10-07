#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Align reads with BWA MEM
process bwa {
    container "${params.container__bwa}"
    tag "cpu_limited"
    
    input:
    tuple val(specimen), path(R1), path(R2)
    path ref

    output:
    tuple val(specimen), path("aligned.bam"), emit: bam

    script:
    template 'bwa.sh'

}

workflow align_wf{

    take:
    reads_ch
    ref

    main:
    bwa(reads_ch, ref)

    emit:
    bam = bwa.out.bam

}