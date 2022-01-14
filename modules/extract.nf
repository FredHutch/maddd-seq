#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

include { merge_bam } from './bwa' addParams(
    subfolder: '6_filtered_SSC/reads',
    file_label: 'adduct.reads',
    publish: true
)

// Extract reads from families which contain adducts
process adduct_reads {
    container "${params.container__bwa}"
    label "mem_medium"
    
    input:
    tuple val(specimen), val(shard_ix), path(bam), path(read_families), path(adduct_families)

    output:
    tuple val(specimen), path("adduct.reads.bam")

    script:
    template 'adduct_reads.sh'

}


workflow extract_wf {

    take:
    bam_shard_ch
    // tuple val(specimen), val(shard_ix), path(bam)
    family_shard_ch
    // tuple val(specimen), val(shard_ix), path("families.csv.gz")
    adduct_family_ch
    // tuple val(specimen), path(adduct_families)

    main:

    adduct_family_ch
        .cross(
            bam_shard_ch
                .join(
                    family_shard_ch,
                    by: [0, 1]
                )
        )
        .map {
            it -> [
                it[0][0], // specimen
                it[1][1], // shard_ix
                it[1][2], // aligned.bam
                it[1][3], // families.csv.gz
                it[0][1]  // adduct.families.txt.gz
            ]
        } | adduct_reads

    // Join the BAM files across all shards
    merge_bam(
        adduct_reads.out.groupTuple()
    )
}