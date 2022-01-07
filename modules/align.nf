#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Import the BWA alignment process with two distinct aliases
include { bwa as align_bwa } from './bwa'
include { bwa as realign_bwa } from './bwa'

// Import the process used to extract the positions of each read
include { extract_positions } from './family'

// Break up the unaligned reads for each specimen into shards for processing
process shard {
    container "${params.container__pandas}"
    label "io_limited"
    
    input:
    tuple val(specimen), path(R1), path(R2), path("barcode_corrections.csv.gz")

    output:
    tuple val(specimen), path("shard.*.R1.fastq.gz"), path("shard.*.R2.fastq.gz")

    script:
    template 'shard.py'

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

// Trim the reads so that they do not overhang the end of the fragment
process trim_overhang {
    container "${params.container__pandas}"
    label "cpu_medium"
    
    input:
    tuple val(specimen), val(shard_ix), path("untrimmed.bam"), path("read_positions.csv.gz")

    output:
    tuple val(specimen), val(shard_ix), path("${specimen}_${shard_ix}_R1.fastq.gz"), path("${specimen}_${shard_ix}_R2.fastq.gz")

    script:
    template 'trim_overhang.sh'

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
    barcodes_csv_ch
    target_regions_bed_path

    main:

    // Break up the unaligned reads into shards which
    // each contain complete sets of barcodes
    shard(
        reads_ch.join(
            barcodes_csv_ch
        )
    )
    // output:
    // tuple val(specimen), path("shard.*.R1.fastq.gz"), path("shard.*.R2.fastq.gz")

    // The output of shard() needs to be transformed to
    // tuple val(specimen), val(shard_ix), path(R1), path(R2)
    shard_ch = shard.out.transpose().map {
        [it[0], it[1].name.replaceAll('.R1.fastq.gz', ''), it[1], it[2]]
    }

    // Align all of the reads
    align_bwa(shard_ch, ref)

    // If the user specified a --target_regions_bed
    if ( target_regions_bed_path ){

        // Get that file
        target_regions_bed = file(target_regions_bed_path)

        // Only keep alignments which overlap this region
        filter_target_regions(
            align_bwa.out.bam,
            target_regions_bed
        )

        bam_ch = filter_target_regions.out.bam
    } else{
        bam_ch = align_bwa.out.bam
    }

    // Extract the positions of each aligned read to enable the trim_overhang method below
    extract_positions(
        bam_ch
    )

    // Merge together the position information CSV with the BAM, using the first two
    // variables in each tuple to join (specimen and shard)

    // The code below is a little ugly. What it's doing is taking the first two variables
    // in each typle and packing them in a nested tuple. This transformation is performed
    // on both of the channels (bam_ch and extract_positions.out), and then the reverse
    // transformation is performed on the resulting channel to give it the expected structure
    // going into trim_overhang.
    bam_ch.map {
        [[it[0], it[1]], it[2]]
    }.join(
        extract_positions.out.map {
            [[it[0], it[1]], it[2]]
        }
    ).map {
        [it[0][0], it[0][1], it[1], it[2]]
    }.set {
        bam_positions_ch
    }

    // Trim the reads so that they do not overhang the end of the fragment
    trim_overhang(
        bam_positions_ch
    )

    // Realign the trimmed reads to the genome
    realign_bwa(trim_overhang.out, ref)

    // Count up the number of aligned reads to each contig per shard
    flagstats(realign_bwa.out)

    // Join flagstats across shards, for each specimen
    join_flagstats(flagstats.out.groupTuple())

    // Assemble a report on the number of reads mapped per specimen
    multiqc_flagstats(join_flagstats.out.toSortedList())

    emit:
    bam = realign_bwa.out

}