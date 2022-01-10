#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Import the BWA alignment process with two distinct aliases
include { bwa as align_bwa } from './bwa'
include { bwa as realign_bwa } from './bwa'

// Import the process used to extract the positions of each read
include { extract_positions } from './family'

// Break up the corrected barcodes for each specimen into shards
process shard_barcodes {
    container "${params.container__pandas}"
    label "io_limited"
    
    input:
    tuple val(specimen), path("barcode_corrections.csv.gz")

    output:
    tuple val(specimen), path("shard.*.barcode_corrections.csv.gz")

    script:
    template 'shard_barcodes.py'

}

// Break up the unaligned reads for each specimen into shards for processing
process shard_reads {
    container "${params.container__pandas}"
    label "io_limited"
    
    input:
    tuple val(specimen), val(shard_ix), path(R1), path(R2), path(shard_barcodes_csv)

    output:
    tuple val(specimen), val(shard_ix), path("${shard_ix}.R1.fastq.gz"), path("${shard_ix}.R2.fastq.gz")

    script:
    template 'shard_reads.py'

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

    // Break up the barcodes into shards which can be analyzed in parallel
    shard_barcodes(
        barcodes_csv_ch
    )

    // Break up the unaligned reads for each specimen into those shards so
    // that each file contains the complete set for each barcode
    shard_reads(
        reads_ch
            .join(shard_barcodes.out)
            .transpose()
            .map {
                [
                    it[0],                   // Specimen ID
                    it[3].name.replaceAll(   // Shard index
                        '.barcode_corrections.csv.gz',
                        ''
                    ),
                    it[1],                   // R1 FASTQ
                    it[2],                   // R2 FASTQ
                    it[3]                    // Barcode list
                ]
            }
    )
    // output:
    // tuple val(specimen), val(shard_ix), path("${shard_ix}.R1.fastq.gz"), path("${shard_ix}.R2.fastq.gz")

    // Align all of the reads
    align_bwa(shard_reads.out, ref)

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