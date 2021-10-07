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

// Count up the number of aligned reads
process flagstats {
    container "${params.container__bwa}"
    publishDir "${params.output}/4_aligned/${specimen}/", mode: 'copy', overwrite: true
    
    input:
    tuple val(specimen), path(bam)

    output:
    file "${specimen}.flagstats"

    script:
    template 'flagstats.sh'

}

// Combine all flagstats data into a single report
process multiqc_flagstats {
    container "${params.container__multiqc}"
    publishDir "${params.output}/4_aligned/", mode: 'copy', overwrite: true
    
    input:
    path "*"

    output:
    path "multiqc_report.html"

    script:
    template 'multiqc.sh'

}

workflow align_wf{

    take:
    reads_ch
    ref

    main:

    // Align all of the reads
    bwa(reads_ch, ref)

    // Count up the number of aligned reads to each contig
    flagstats(bwa.out.bam)

    // Assemble a report on the number of reads mapped per specimen
    multiqc_flagstats(flagstats.out.toSortedList())

    emit:
    bam = bwa.out.bam

}