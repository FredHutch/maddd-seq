#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Break up the unaligned reads for each specimen into shards for processing
process shard {
    container "${params.container__pandas}"
    label "io_limited"
    
    input:
    tuple val(specimen), path(R1), path(R2)

    output:
    tuple val(specimen), path("shard.*.R1.fastq.gz"), path("shard.*.R2.fastq.gz")

    script:
    template 'shard.py'

}

// Align reads with BWA MEM
process bwa {
    container "${params.container__bwa}"
    label "cpu_medium"
    
    input:
    tuple val(specimen), val(shard_ix), path(R1), path(R2)
    path ref

    output:
    tuple val(specimen), val(shard_ix), path("aligned.bam"), emit: bam

    script:
    template 'bwa.sh'

}

// Filter all alignments to those which overlap a target region
process filter_target_regions {
    container "${params.container__bwa}"
    label "cpu_medium"
    
    input:
    tuple val(specimen), val(shard_ix), path("unmasked.bam")
    path target_regions_bed

    output:
    tuple val(specimen), val(shard_ix), path("aligned.bam"), emit: bam

    script:
    template 'filter_target_regions.sh'

}

// Count up the number of aligned reads
process flagstats {
    container "${params.container__bwa}"
    label "io_limited"
    
    input:
    tuple val(specimen), val(shard_ix), path(bam)

    output:
    tuple val(specimen), path("${shard_ix}.flagstats")

    script:
    template 'flagstats.sh'

}

// Join flagstats across shards per specimen
process join_flagstats {
    container "${params.container__bwa}"
    publishDir "${params.output}/4_aligned/${specimen}/", mode: 'copy', overwrite: true
    label "io_limited"
    
    input:
    tuple val(specimen), path("*")  // "${shard_ix}.flagstats"

    output:
    path "${specimen}.flagstats"

    script:
    template 'join_flagstats.py'

}

// Combine all flagstats data into a single report
process multiqc_flagstats {
    container "${params.container__multiqc}"
    publishDir "${params.output}/4_aligned/", mode: 'copy', overwrite: true
    label "io_limited"
    
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

    // Break up the unaligned reads into shards which
    // each contain complete sets of barcodes
    shard(
        reads_ch
    )
    // output:
    // tuple val(specimen), path("shard.*.R1.fastq.gz"), path("shard.*.R2.fastq.gz")

    // The output of shard() needs to be transformed to
    // tuple val(specimen), val(shard_ix), path(R1), path(R2)
    shard_ch = shard.out.transpose().map {
        [it[0], it[1].name.replaceAll('.R1.fastq.gz', ''), it[1], it[2]]
    }

    // Align all of the reads
    bwa(shard_ch, ref)

    // If the user specified a --target_regions_bed
    if ( params.target_regions_bed ){

        // Get that file
        target_regions_bed = file(params.target_regions_bed)

        // Only keep alignments which overlap this region
        filter_target_regions(
            bwa.out.bam,
            target_regions_bed
        )

        bam_ch = filter_target_regions.out.bam
    } else{
        bam_ch = bwa.out.bam
    }

    // Count up the number of aligned reads to each contig per shard
    flagstats(bwa.out.bam)

    // Join flagstats across shards, for each specimen
    join_flagstats(flagstats.out.groupTuple())

    // Assemble a report on the number of reads mapped per specimen
    multiqc_flagstats(join_flagstats.out.toSortedList())

    emit:
    bam = bwa.out.bam

}