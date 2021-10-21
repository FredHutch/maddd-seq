#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Filter SSCs based on the depth of sequencing
process filter_ssc_depth {
    container "${params.container__pandas}"
    publishDir "${params.output}/7_filtered_SSC/${specimen}/alignments/", mode: 'copy', overwrite: true, pattern: "*.bam"
    
    input:
    tuple val(specimen), path("unfiltered.POS.SSC.bam"), path("unfiltered.NEG.SSC.bam"), path("unfiltered.SSC.details.csv.gz")

    output:
    tuple val(specimen), path("POS.SSC.bam"), path("NEG.SSC.bam"), path("SSC.details.csv.gz"), optional: true

    script:
    template 'filter_ssc_depth.py'

}


// Parse the SSC data
process parse_ssc {
    container "${params.container__pandas}"
    publishDir "${params.output}/7_filtered_SSC/${specimen}/stats/", mode: 'copy', overwrite: true, pattern: "*csv.gz"
    
    input:
    tuple val(specimen), path("POS.SSC.bam"), path("NEG.SSC.bam"), path("SSC.details.csv.gz")

    output:
    tuple val(specimen), path("total.json.gz"), emit: json
    path "summary.json.gz"
    path "*.csv.gz", emit: csv

    script:
    template 'parse_ssc.py'

}

// Format the DSC data as BAM
process format_dsc {
    container "${params.container__pandas}"
    publishDir "${params.output}/7_filtered_SSC/${specimen}/alignments/", mode: 'copy', overwrite: true
    
    input:
    tuple val(specimen), path("POS.SSC.bam"), path("NEG.SSC.bam"), path("SSC.details.csv.gz")

    output:
    tuple val(specimen), path("DSC.bam")

    script:
    template 'format_dsc.py'

}

// Format the output as VCF
process format_vcf {
    container "${params.container__bcftools}"
    publishDir "${params.output}/7_filtered_SSC/${specimen}/alignments/", mode: 'copy', overwrite: true
    
    input:
    tuple val(specimen), path("DSC.bam")
    path ref

    output:
    file "DSC.vcf.gz"

    script:
    template 'format_vcf.sh'

}

// Format details about all SSCs as CSV
process format_ssc_csv {
    container "${params.container__pandas}"
    publishDir "${params.output}/7_filtered_SSC/${specimen}/stats/", mode: 'copy', overwrite: true
    
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
    publishDir "${params.output}/7_filtered_SSC/plots/", mode: 'copy', overwrite: true
    
    input:
    path "*"

    output:
    file "*.pdf"

    script:
    template 'make_plots.py'

}

workflow variants_wf{

    take:
    bam_ch
    // tuple val(specimen), path("POS.SSC.bam"), path("NEG.SSC.bam"), path("SSC.details.csv.gz")
    genome_ref

    main:

    // Filter the SSC data based on --min_reads
    filter_ssc_depth(
        bam_ch
    )
    // output:
    // tuple val(specimen), path("POS.SSC.bam"), path("NEG.SSC.bam"), path("SSC.details.csv.gz")

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

    // Format the output as VCF
    format_vcf(
        format_dsc.out,
        genome_ref
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
            .toSortedList()
    )

}