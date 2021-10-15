#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Break up the unaligned reads for each specimen into shards for processing
process shard {
    container "${params.container__pandas}"
    
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
    tag "cpu_limited"
    
    input:
    tuple val(specimen), val(shard_ix), path(R1), path(R2)
    path ref

    output:
    tuple val(specimen), val(shard_ix), path("aligned.bam"), emit: bam

    script:
    template 'bwa.sh'

}

// Count up the number of aligned reads
process flagstats {
    container "${params.container__bwa}"
    
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
    // tuple val(specimen), val(shard_ix), path(bam)
    shard_ch = shard.out.transpose().map {
        [it[0], it[1].name.replaceAll('.R1.fastq.gz', ''), it[1], it[2]]
    }

    // Align all of the reads
    bwa(shard_ch, ref)

    // Count up the number of aligned reads to each contig per shard
    flagstats(bwa.out.bam)

    // Join flagstats across shards, for each specimen
    join_flagstats(flagstats.out.groupTuple())

    // Assemble a report on the number of reads mapped per specimen
    multiqc_flagstats(join_flagstats.out.toSortedList())

    emit:
    bam = bwa.out.bam

}