#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Import the BWA alignment process with two distinct aliases
include { bwa as align_bwa } from './bwa'
include { bwa as realign_bwa } from './bwa'
include { merge_bam } from './bwa' addParams(
    subfolder: '4_aligned/reads',
    file_label: 'aligned',
    publish: true
)

// Import the process used to extract the positions of each read
include { extract_positions } from './family'

// Break up the corrected barcodes for each specimen into shards
process shard_barcodes {
    container "${params.container__pandas}"
    publishDir "${params.output}/4_aligned/shard_intermediate/", mode: 'copy', enabled: params.save_intermediates
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

// Merge paired FASTQ files across shards
process trim_overhang_join_shards {
    container "${params.container__bwa}"
    publishDir "${params.output}/4_aligned/trim_overhang/${specimen}/", mode: 'copy', overwrite: true
    label "io_limited"
    
    input:
    tuple val(specimen), path("R1/"), path("R2/")

    output:
    tuple path("${specimen}_R1.fastq.gz"), path("${specimen}_R2.fastq.gz")

    script:
    template 'join_shards_fastq.sh'

}

// Count up the number of aligned reads
process flagstats {
    container "${params.container__bwa}"
    publishDir "${params.output}/4_aligned/flagstats_intermediate/${shard_ix}/", mode: 'copy', enabled: params.save_intermediates
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
    publishDir "${params.output}/4_aligned/join_flagstats_intermediate/", mode: 'copy', enabled: params.save_intermediates
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

    // Extract the positions of each aligned read to enable the trim_overhang method below
    extract_positions(
        align_bwa.out.bam
    )

    // Merge together the position information CSV with the BAM, using the first two
    // variables in each tuple to join (specimen and shard)

    // The code below is a little ugly. What it's doing is taking the first two variables
    // in each typle and packing them in a nested tuple. This transformation is performed
    // on both of the channels (align_bwa.out.bam and extract_positions.out), and then the reverse
    // transformation is performed on the resulting channel to give it the expected structure
    // going into trim_overhang.
    align_bwa.out.bam.map {
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

    // Publish the overhang-trimmed reads
    trim_overhang_join_shards(
        trim_overhang
            .out
            .map {
                it -> [it[0], it[2], it[3]]
            }.groupTuple()
    )

    // Realign the trimmed reads to the genome
    realign_bwa(trim_overhang.out, ref)

    // Publish the aligned reads after merging the BAMs from each shard
    merge_bam(
        realign_bwa
            .out
            .bam
            .map {
                it -> [it[0], it[2]]
            }.groupTuple()
    )

    // Count up the number of aligned reads to each contig per shard
    flagstats(realign_bwa.out)

    // Join flagstats across shards, for each specimen
    join_flagstats(flagstats.out.groupTuple())

    // Assemble a report on the number of reads mapped per specimen
    multiqc_flagstats(join_flagstats.out.toSortedList())

    emit:
    bam = realign_bwa.out

}