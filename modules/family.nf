#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Break up the aligned reads for each specimen into shards for processing
process shard {
    container "${params.container__pysam}"
    
    input:
    tuple val(specimen), path(bam), path(barcodes_csv_gz)

    output:
    tuple val(specimen), path("shard.*.bam")

    script:
    template 'shard.py'

}

// Extract the ID, chromosome, position, orientation, and R1/R2 for each alignment
process extract_positions {
    container "${params.container__pysam}"
    
    input:
    tuple val(specimen), val(shard_ix), path(bam)

    output:
    tuple val(specimen), val(shard_ix), path("read_positions.csv.gz")

    script:
    template 'extract_positions.py'

}

workflow family_wf{

    take:
    input_ch
    // tuple val(specimen), path(bam), path(barcodes_csv_gz)

    main:

    // Break up the aligned BAM into shards which
    // each contain a set of barcodes
    shard(
        input_ch
    )
    // output:
    // tuple val(specimen), path("shard.*.bam")

    // The output of shard() needs to be transformed to
    // tuple val(specimen), val(shard_ix), path(bam)
    shard_ch = shard.out.transpose().map {
        [it[0], it[1].name.replaceAll('.bam', ''), it[1]]
    }

    // Extract the ID, position, orientation, chromosome,
    // and R1/R2 for each alignment
    extract_positions(
        shard_ch
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