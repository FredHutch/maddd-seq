#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Filter SSCs based on the depth of sequencing
process filter_ssc_depth {
    container "${params.container__pandas}"
    publishDir "${params.output}/6_filtered_SSC/${specimen}/filter_ssc_depth_intermediate/", mode: 'copy', enabled: params.save_intermediates
    label "io_limited"
    
    input:
    tuple val(specimen), path("unfiltered.POS.SSC.bam"), path("unfiltered.NEG.SSC.bam"), path("unfiltered.SSC.details.csv.gz")

    output:
    tuple val(specimen), path("POS.SSC.bam"), path("NEG.SSC.bam"), path("SSC.details.csv.gz"), optional: true

    script:
    template 'filter_ssc_depth.py'

}

// Sort and index SSC BAM files
process index_ssc {
    container "${params.container__bwa}"
    publishDir "${params.output}/6_filtered_SSC/${specimen}/alignments/", mode: 'copy', overwrite: true
    label "io_limited"
    
    input:
    tuple val(specimen), val(output_filename), path("unfiltered.bam")

    output:
    tuple val(specimen), path("${output_filename}"), path("${output_filename}.bai")

    script:
    template 'index_ssc.sh'

}


// Parse the SSC data
process parse_ssc {
    container "${params.container__pandas}"
    publishDir "${params.output}/6_filtered_SSC/${specimen}/stats/", mode: 'copy', overwrite: true, pattern: "*csv.gz"
    publishDir "${params.output}/6_filtered_SSC/${specimen}/stats/", mode: 'copy', overwrite: true, pattern: "*adduct.gtf"
    label "io_limited"
    
    input:
    tuple val(specimen), path("POS.SSC.bam"), path("NEG.SSC.bam"), path("SSC.details.csv.gz")

    output:
    tuple val(specimen), path("total.json.gz"), emit: json
    path "summary.json.gz"
    path "*.csv.gz", emit: csv
    tuple val(specimen), path("${specimen}.adduct.families.txt.gz"), emit: adduct_families
    path "${specimen}.adduct.gtf"

    script:
    template 'parse_ssc.sh'

}

// Format the DSC data as BAM
process format_dsc {
    container "${params.container__pandas}"
    publishDir "${params.output}/6_filtered_SSC/${specimen}/format_dsc_intermediate/", mode: 'copy', enabled: params.save_intermediates
    label "io_limited"
    
    input:
    tuple val(specimen), path("POS.SSC.bam"), path("NEG.SSC.bam"), path("SSC.details.csv.gz")

    output:
    tuple val(specimen), path("DSC.bam")

    script:
    template 'format_dsc.py'

}

// Sort and index the BAM file for each DSC
process index_dsc {
    container "${params.container__bwa }"
    publishDir "${params.output}/6_filtered_SSC/${specimen}/alignments/", mode: 'copy', overwrite: true
    label "io_limited"
    
    input:
    tuple val(specimen), path("unsorted.DSC.bam")

    output:
    tuple val(specimen), path("DSC.bam")
    path "DSC.bam.bai"

    script:
    template 'index_dsc.sh'

}

// Format the align DSC reads as VCF
process format_vcf {
    container "${params.container__bcftools}"
    publishDir "${params.output}/6_filtered_SSC/${specimen}/alignments/", mode: 'copy', overwrite: true
    label "io_limited"
    
    input:
    tuple val(specimen), path("DSC.bam")
    path ref

    output:
    file "DSC.vcf.gz"

    script:
    template 'format_vcf.sh'

}

// Format the aligned DSC reads as pileup
process format_pileup {
    container "${params.container__bwa}"
    publishDir "${params.output}/6_filtered_SSC/${specimen}/alignments/", mode: 'copy', overwrite: true
    label "io_limited"
    
    input:
    tuple val(specimen), path("DSC.bam")
    path ref

    output:
    tuple val(specimen), path("DSC.pileup.gz")

    script:
    template 'format_pileup.sh'

}

// Format the align DSC reads as TSV
process format_tsv {
    container "${params.container__pandas}"
    publishDir "${params.output}/6_filtered_SSC/${specimen}/stats/", mode: 'copy', overwrite: true
    label "io_limited"
    
    input:
    tuple val(specimen), path(pileup)

    output:
    file "${specimen}.DSC.tsv.gz"

    script:
    """#!/bin/bash
    format_tsv.py $pileup ${specimen}.DSC.tsv.gz
    """

}

// Format details about all SSCs as CSV
process format_ssc_csv {
    container "${params.container__pandas}"
    publishDir "${params.output}/6_filtered_SSC/${specimen}/stats/", mode: 'copy', overwrite: true
    label "io_limited"
    
    input:
    tuple val(specimen), path("total.json.gz"), path("SSC.details.csv.gz")

    output:
    file "${specimen}.SSC.csv.gz"

    script:
    template 'format_ssc_csv.py'

}

// Make plots
process make_plots {
    container "${params.container__python_plotting}"
    publishDir "${params.output}/6_filtered_SSC/plots/", mode: 'copy', overwrite: true, pattern: "*.pdf"
    publishDir "${params.output}/6_filtered_SSC/tables/", mode: 'copy', overwrite: true, pattern: "*.csv"
    label "io_limited"
    
    input:
    path "*"

    output: 
    file "*"

    script:
    template 'make_plots.sh'

}

workflow variants_wf{

    take:
    bam_ch
    // tuple val(specimen), path("POS.SSC.bam"), path("NEG.SSC.bam"), path("SSC.details.csv.gz")
    genome_ref
    // CSV with the number of reads per barcode, for each specimen
    barcode_counts

    main:

    // Filter the SSC data based on --min_reads
    filter_ssc_depth(
        bam_ch
    )
    // output:
    // tuple val(specimen), path("POS.SSC.bam"), path("NEG.SSC.bam"), path("SSC.details.csv.gz")

    // To publish the SSC bam files for output, they must be sorted and indexed
    // Note in this case that the files which are published do not have the exact order
    // as the files which are used to format the DSC data, below. This should not impact
    // the results, but it may be important for troubleshooting purposes
    // NOTE: the index_ssc process will be executed on both the POS and NEG SSC data independently
    index_ssc(
        filter_ssc_depth
            .out
            .map {[it[0], "POS.SSC.bam", it[1]]}
            .mix(
            filter_ssc_depth
                .out
                .map {[it[0], "NEG.SSC.bam", it[2]]}
            )
    )

    // Parse the SSC data in order to:
    //   - Construct DSCs BAMs which combine both strands
    //   - Call adducts from the SSC data
    //   - Call SNPs from the SSC data
    //   - Summarize the total number of adducts and SNPs
    //   - Summarize the number of adducts and SNPs per chromosome
    //   - Summarize the number of adducts and SNPs per position within each read
    parse_ssc(
        filter_ssc_depth.out
    )

    // Format the DSC data as BAM
    format_dsc(
        filter_ssc_depth.out
    )

    // Sort the DSC BAM file
    index_dsc(
        format_dsc.out
    )

    // Format the output as VCF
    format_vcf(
        index_dsc.out[0],
        genome_ref
    )

    // Make a pileup file
    format_pileup(
        index_dsc.out[0],
        genome_ref
    )

    // Format the output as TSV
    format_tsv(
        format_pileup.out
    )

    // Format details about all SSCs as CSV
    format_ssc_csv(
        // Format the input to this as the JSON with mutations per 
        // family, as well as the number of reads per SSC
        parse_ssc
            .out
            .json
            .join(
                filter_ssc_depth
                    .out
                    .map {
                        [it[0], it[3]]
                    }
            )
    )

    // Make a series of plots across all specimens
    make_plots(
        format_ssc_csv
            .out
            .mix(
                parse_ssc
                    .out
                    .csv
                    .flatten()
            )
            .mix(
                barcode_counts
            )
            .mix(
                bam_ch.map {
                    it -> it[3]
                }
            )
            .toSortedList()
    )

    emit:
    adduct_families = parse_ssc.out.adduct_families

}