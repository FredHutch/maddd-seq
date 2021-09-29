#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

workflow align_wf{

    take:
    bam_ch

    main:
    align(bam_ch)

    emit:
    aligned_bam = align.out

}