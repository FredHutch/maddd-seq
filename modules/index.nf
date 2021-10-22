#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2


// Index a genome with BWA MEM
process bwa_index {
    container "${params.container__bwa}"
    tag "cpu_limited"
    
    input:
    path genome_fasta

    output:
    // Output all files which were created with this command
    path "*"

    script:
    template 'bwa_index.sh'

}

// Run RepeatMasker on a genome FASTA
process repeatmasker {
    container "${params.container__repeatmasker}"
    tag "cpu_limited"
    
    input:
    path genome_fasta

    output:
    // Output all files which were created with this command
    path "*"

    script:
    template 'repeatmasker.sh'

}
