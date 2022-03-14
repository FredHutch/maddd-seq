#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2


// Extract reads from families which contain adducts
process adduct_reads {
    container "${params.container__bwa}"
    label "mem_medium"
    
    input:
    tuple val(specimen), val(filtering), path(adduct_families), val(shard_ix), path(bam), path(read_families)

    output:
    tuple val(specimen), val(filtering), path("adduct.reads.bam")

    script:
    template 'adduct_reads.sh'

}


// Merge a collection of alignments in sorted BAM format
process merge_adduct_reads_bam {
    container "${params.container__bwa}"
    publishDir "${params.output}/6_filtered_SSC/${specimen}/${filtering}/", mode: 'copy', overwrite: true
    label "cpu_medium"
    
    input:
    tuple val(specimen), val(filtering), path("input_bam/*.bam")

    output:
    path "*"

    script:
    template 'merge_adduct_reads_bam.sh'

}


workflow extract_wf {

    take:
    bam_shard_ch
    // tuple val(specimen), val(shard_ix), path(bam)
    family_shard_ch
    // tuple val(specimen), val(shard_ix), path("families.csv.gz")
    adduct_family_ch
    // tuple val(specimen), val(filtering), path(adduct_families)

    main:

    adduct_family_ch
        .combine(
            bam_shard_ch
                .join(
                    family_shard_ch,
                    by: [0, 1]
                ),
            by: 0
        ) | adduct_reads

    // Join the BAM files across all shards
    merge_adduct_reads_bam(
        adduct_reads.out.groupTuple(
            by: [0, 1]
        )
    )
}