#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Extract the ID, chromosome, position, orientation, and R1/R2 for each alignment
process extract_positions {
    container "${params.container__pandas}"
    publishDir "${params.output}/5_all_SSC/${specimen}/extract_positions_intermediate/", mode: 'copy', enabled: params.save_intermediates
    label "io_limited"
    
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
    publishDir "${params.output}/5_all_SSC/${specimen}/assign_families_intermediate/", mode: 'copy', enabled: params.save_intermediates
    label "io_limited"
    
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
    publishDir "${params.output}/5_all_SSC/${specimen}/make_ssc_intermediate/", mode: 'copy', enabled: params.save_intermediates
    label "io_limited"
    
    input:
    tuple val(specimen), val(shard_ix), path("families.csv.gz"), path(bam)

    output:
    tuple val(specimen), val(shard_ix), path("${shard_ix}.FWD.R1.fastq.gz"), path("${shard_ix}.FWD.R2.fastq.gz"), path("${shard_ix}.REV.R1.fastq.gz"), path("${shard_ix}.REV.R2.fastq.gz"), path("${shard_ix}.stats.csv.gz"), optional: true

    script:
    template 'make_ssc.py'

}

// Re-align the SSC sequences against the reference
process align_ssc {
    container "${params.container__bwa}"
    publishDir "${params.output}/5_all_SSC/${specimen}/align_ssc_intermediate/", mode: 'copy', enabled: params.save_intermediates
    label "cpu_medium"
    
    input:
    tuple val(specimen), path("SSC_FWD_R1/*"), path("SSC_FWD_R2/*"), path("SSC_REV_R1/*"), path("SSC_REV_R2/*")
    path ref

    output:
    tuple val(specimen), path("POS.SSC.bam"), path("NEG.SSC.bam")

    script:
    template 'align_ssc.sh'

}

// Filter out any re-aligned sequences whose anchor (outermost) coordinates
// have changed after the realignment of SSCs
process filter_ssc_position {
    container "${params.container__pandas}"
    publishDir "${params.output}/5_all_SSC/${specimen}/", mode: 'copy', overwrite: true, pattern: "*.csv.gz"
    label "mem_verylarge"
    
    input:
    tuple val(specimen), path("realigned.POS.SSC.bam"), path("realigned.NEG.SSC.bam"), path("SSC_STATS/*.csv.gz")

    output:
    tuple val(specimen), path("POS.SSC.bam"), path("NEG.SSC.bam"), path("${specimen}.unfiltered.SSC.details.csv.gz")

    script:
    template 'filter_ssc_position.py'

}

// Sort and index the SSC BAM files
process sort_ssc_bam {
    container "${params.container__bwa}"
    publishDir "${params.output}/5_all_SSC/${specimen}/", mode: 'copy', overwrite: true, pattern: "*.SSC.ba*"
    label "mem_medium"
    
    input:
    tuple val(specimen), path("unsorted.POS.SSC.bam"), path("unsorted.NEG.SSC.bam"), path("SSC.details.csv.gz")

    output:
    file "*"

    script:
    template 'sort_ssc_bam.sh'

}

workflow family_wf{

    take:
    bam_shard_ch
    // tuple val(specimen), val(shard_ix), path(bam)
    ref
    // Pre-compiled genome reference

    main:

    // Extract the ID, position, orientation, chromosome,
    // and R1/R2 for each alignment
    extract_positions(
        bam_shard_ch
    )

    // Group reads into families which share the same
    // barcode and alignment position. The family ID
    // will be encoded in the attached CSV by this step
    assign_families(
        extract_positions.out
    )
    // input:
    //   tuple val(specimen), val(shard_ix), path("read_positions.csv.gz")
    // output:
    //   tuple val(specimen), val(shard_ix), path("families.csv.gz")

    // Compute the SSC sequences at the FASTQ level
    // Note that this requires reads from *both* strands for each family
    make_ssc(
        assign_families
            .out
            .combine(bam_shard_ch, by: [0, 1])
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
    //   tuple val(specimen), path(POS_BAM), path(NEG_BAM)

    // Filter out any re-aligned sequences whose anchor (outermost) coordinates
    // have changed after the realignment of SSCs
    SSC_STATS = make_ssc.out.map { [it[0], it[6]] }.groupTuple()
    filter_ssc_position(
        align_ssc
            .out
            .join(SSC_STATS)
    )
    // input:
    //   tuple val(specimen), path(POS_BAM), path(NEG_BAM), path(SSC_STATS/*csv.gz)
    // output:
    //   tuple val(specimen), path("POS.SSC.bam"), path("NEG.SSC.bam"), path("SSC.details.csv.gz")

    // Sort and index the SSC BAM files
    sort_ssc_bam(
        filter_ssc_position.out
    )

    emit:
    bam = filter_ssc_position.out
    families = assign_families.out

}