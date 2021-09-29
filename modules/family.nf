#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

workflow family_wf{

    take:
    aligned_bam_ch

    main:

    // Break up the aligned BAM into shards which
    // each contain a set of barcodes
    shard(
        aligned_bam_ch
    )

    // Group reads into families which share the same
    // barcode and alignment position. The family ID
    // will be encoded in the BAM header by this step
    assign_families(shard.out)

    // Compute the SSC sequences
    make_ssc(assign_families.out)
    //output: tuple val(specimen), path(ssc_bam)

    // Compute the DSC sequences
    make_dsc(make_ssc.out)
    //output: tuple val(specimen), path(dsc_bam)

    // Summarize the characteristics of each SSC
    summarize_ssc(make_ssc.out)

    // Summarize the characteristics of each DSC
    summarize_dsc(make_dsc.out)

    emit:
    // tuple val(specimen), path(ssc_bam), path(dsc_bam)
    bam = make_ssc.out.join(
        make_dsc.out,
        by: [0, 1]
    )

    // tuple val(specimen), path(ssc_summary_csv)
    ssc_summary = summarize_ssc.out

    // tuple val(specimen), path(dsc_summary_csv)
    dsc_summary = summarize_dsc.out

}