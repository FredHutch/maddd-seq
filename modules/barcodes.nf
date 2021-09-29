#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

workflow barcodes_wf{

    take:
    fastq_ch

    main:

    // Remove the barcodes from each end of the read
    // and move to the header of an unaligned BAM
    clip_barcodes(fastq_ch)

    // Count up the number of reads with each barcode
    count_barcodes(clip_barcodes.out)

    // Perform error correction on the barcodes
    correct_barcode_errors(count_barcodes.out)

    // Modify the barcode in the unaligned BAM to
    // the corrected sequence
    update_barcodes(
        fastq_ch.combine(
            correct_barcode_errors.out
        )
    )

    emit:
    bam = update_barcodes.out

}