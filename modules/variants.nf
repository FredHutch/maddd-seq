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
    publishDir "${params.output}/6_filtered_SSC/${specimen}/", mode: 'copy', overwrite: true
    label "mem_medium"
    
    input:
    tuple val(specimen), path("POS.SSC.bam"), path("NEG.SSC.bam"), path("SSC.details.csv.gz")

    output:
    path "*"
    tuple val(specimen), path("*/*.DSC.bam"), emit: dsc_bam
    tuple val(specimen), path("*/*.SSC.POS.bam"), emit: ssc_pos_bam
    tuple val(specimen), path("*/*.SSC.NEG.bam"), emit: ssc_neg_bam
    tuple val(specimen), path("*/*.csv"), emit: csv
    tuple val(specimen), path("*/*.json.gz"), emit: json
    tuple val(specimen), path("*/*.adduct.families.txt.gz"), emit: adduct_families

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
    publishDir "${params.output}/6_filtered_SSC/${specimen}/${filtering}/", mode: 'copy', overwrite: true
    label "io_limited"
    
    input:
    tuple val(specimen), val(filtering), path(bam)
    path ref

    output:
    file "${bam.name.replaceAll(".bam", "")}.vcf.gz"

    script:
    template 'format_vcf.sh'

}

// Format the aligned DSC reads as pileup
process format_pileup {
    container "${params.container__bwa}"
    publishDir "${params.output}/6_filtered_SSC/${specimen}/${filtering}/", mode: 'copy', overwrite: true
    label "io_limited"
    
    input:
    tuple val(specimen), val(filtering), path(bam)
    path ref

    output:
    tuple val(specimen), val(filtering), path("${bam.name.replaceAll(".bam", "")}.pileup.gz")

    script:
    template 'format_pileup.sh'

}

// Format the align DSC reads as TSV
process format_tsv {
    container "${params.container__pandas}"
    publishDir "${params.output}/6_filtered_SSC/${specimen}/${filtering}/", mode: 'copy', overwrite: true
    label "io_limited"
    
    input:
    tuple val(specimen), val(filtering), path(pileup)

    output:
    file "${pileup.name.replaceAll(".pileup.gz", "")}.tsv.gz"

    script:
    """#!/bin/bash
    format_tsv.py $pileup ${pileup.name.replaceAll(".pileup.gz", "")}.tsv.gz
    """

}

// Format details about all SSCs as CSV
process format_ssc_csv {
    container "${params.container__pandas}"
    publishDir "${params.output}/6_filtered_SSC/${specimen}/${filtering}/", mode: 'copy', overwrite: true
    label "io_limited"
    
    input:
    tuple val(specimen), val(filtering), path("total.json.gz"), path("SSC.details.csv.gz")

    output:
    file "SSC.csv.gz"

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
    //   - Call adducts from the SSC data
    //   - Call SNPs from the SSC data
    //   - Count the number of mutations and adducts per DSC
    //   - Make a table of all mutations and adducts per DSC
    //   - Write out a DSC BAM and SSC BAM for each possible
    //     level of filtering based on the maximum number of
    //     allowable mutations and adducts per DSC
    //   - Make a summary table of the total number of adducts
    //     and SNPs for each of the possible levels of filtering
    // The output of this process will be published to
    //   ${params.output}/6_filtered_SSC/${specimen}/
    // with the suffixes:
    //   ${params.output}/6_filtered_SSC/${specimen}/all/
    //   ${params.output}/6_filtered_SSC/${specimen}/max_variants_{i}/
    // where {i} is a non-negative integer with the maximum
    // number of allowed variants per family (the sum of
    // adducts and variants) in the filtered results
    // Note: All of the files inside each subfolder will have the
    // prefix all or max_variants_{i} so that downstream processes
    // can identify the appropriate subfolder for publishing any
    // additional outputs which will be derived from them
    parse_ssc(
        filter_ssc_depth.out
    )

    // Make a channel which contains each of the BAM files
    // from the parse_ssc output, along with a variable
    // indicating the filtering which was performed on the
    // families based on the maximum number of variants and adducts
    parse_ssc
        .out
        .dsc_bam
        .transpose()
        .map {
            it -> [it[0], it[1].name.replaceAll(".DSC.bam", ""), it[1]]
        }
        .mix(
            parse_ssc
                .out
                .ssc_pos_bam
                .transpose()
                .map {
                    it -> [it[0], it[1].name.replaceAll(".SSC.POS.bam", ""), it[1]]
                }
        )
        .mix(
            parse_ssc
                .out
                .ssc_neg_bam
                .transpose()
                .map {
                    it -> [it[0], it[1].name.replaceAll(".SSC.NEG.bam", ""), it[1]]
                }
        )
        .set { merged_bam_ch }

    // Make a pileup file from each of the BAM files
    format_pileup(
        merged_bam_ch,
        genome_ref
    )

    // Format the pileup data as VCF
    format_vcf(
        merged_bam_ch,
        genome_ref
    )

    // Format the output as TSV from the pileup
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
            .transpose()
            .map {
                it -> [it[0], it[1].name.replaceAll(".json.gz", ""), it[1]]
            }
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