#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Break up the aligned reads for each specimen into shards for processing
process shard {
    container "${params.container__pandas}"
    
    input:
    tuple val(specimen), path(bam), path(barcodes_csv_gz)

    output:
    tuple val(specimen), path("shard.*.bam")

    script:
    template 'shard.py'

}

// Extract the ID, chromosome, position, orientation, and R1/R2 for each alignment
process extract_positions {
    container "${params.container__pandas}"
    
    input:
    tuple val(specimen), val(shard_ix), path(bam)

    output:
    tuple val(specimen), val(shard_ix), path("read_positions.csv.gz")

    script:
    template 'extract_positions.py'

}

// Group reads into families which share the same
// barcode and alignment position. The family ID
// will be encoded in the attached CSV by this step
process assign_families {
    container "${params.container__pandas}"
    
    input:
    tuple val(specimen), val(shard_ix), path("read_positions.csv.gz")

    output:
    tuple val(specimen), val(shard_ix), path("families.csv.gz")

    script:
    template 'assign_families.py'

}

// Compute the SSC sequences at the FASTQ level
process make_ssc {
    container "${params.container__pandas}"
    
    input:
    tuple val(specimen), val(shard_ix), path("families.csv.gz"), path(bam)

    output:
    tuple val(specimen), val(shard_ix), path("FWD.R1.fastq.gz"), path("FWD.R2.fastq.gz"), path("REV.R1.fastq.gz"), path("REV.R2.fastq.gz"), path("${shard_ix}.stats.csv.gz")

    script:
    template 'make_ssc.py'

}

// Re-align the SSC sequences against the reference
process align_ssc {
    container "${params.container__bwa}"
    
    input:
    tuple val(specimen), path("SSC_FWD_R1/*"), path("SSC_FWD_R2/*"), path("SSC_REV_R1/*"), path("SSC_REV_R2/*")
    path ref

    output:
    tuple val(specimen), val(shard_ix), path("POS.SSC.bam"), path("NEG.SSC.bam")

    script:
    template 'align_ssc.sh'

}
workflow family_wf{

    take:
    bam_ch
    // tuple val(specimen), path(bam), path(barcodes_csv_gz)
    ref
    // Pre-compiled genome reference

    main:

    // Break up the aligned BAM into shards which
    // each contain a set of barcodes
    shard(
        bam_ch
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
    // will be encoded in the attached CSV by this step
    assign_families(
        extract_positions.out
    )

    // Compute the SSC sequences at the FASTQ level
    // Note that this requires reads from *both* strands for each family
    make_ssc(
        assign_families
            .out
            .combine(shard_ch, by: [0, 1])
    )
    // input:
    //   tuple val(specimen), val(shard_ix), path(families_csv_gz), path(R1), path(R2)
    // output: 
    //   tuple val(specimen), val(shard_ix), path(SSC_FWD_R1), path(SSC_FWD_R2), path(SSC_REV_R1), path(SSC_REV_R2), path(SSC_STATS_CSV)

    // Join the shards together for each specimen, with each part of the family broken out
    // The complex logic below groups each of the FASTQ files into lists
    FWD_R1 = make_ssc.out.map { [it[0], it[2]] }.groupTuple()
    FWD_R2 = make_ssc.out.map { [it[0], it[3]] }.groupTuple()
    REV_R1 = make_ssc.out.map { [it[0], it[4]] }.groupTuple()
    REV_R2 = make_ssc.out.map { [it[0], it[5]] }.groupTuple()

    // Re-align the SSC consensus sequences against the reference
    // When re-aligning, the FWD_R1 reads will be paired with the REV_R2, and
    // the resulting alignment will be labeled "POS" (as in strand)
    // correspondingly, FWD_R2 reads will be paired with the REV_R1 reads
    // with the label "NEG"
    align_ssc(
        // Filter to the specimen and FASTQ data only, and group by specimen
        FWD_R1
            .join(FWD_R2)
            .join(REV_R1)
            .join(REV_R2),
        // Reference genome
        ref
    )
    // output:
    //   tuple val(specimen), val(shard_ix), path(POS_BAM), path(NEG_BAM)

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