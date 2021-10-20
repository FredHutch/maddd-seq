#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

// Filter SSCs based on the depth of sequencing
process filter_ssc_depth {
    container "${params.container__pandas}"
    publishDir "${params.output}/7_filtered_SSC/${specimen}/", mode: 'copy', overwrite: true
    
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
    publishDir "${params.output}/7_filtered_SSC/${specimen}/", mode: 'copy', overwrite: true
    
    input:
    tuple val(specimen), path("POS.SSC.bam"), path("NEG.SSC.bam"), path("SSC.details.csv.gz")

    output:
    tuple val(specimen), path("total.json.gz"), emit: json
    file "summary.json.gz"
    file "by_chr.csv.gz"
    file "snps_by_base.csv.gz"
    file "adducts_by_base.csv.gz"
    file "by_read_position.csv.gz"

    script:
    template 'parse_ssc.py'

}

// Format the DSC data as BAM
process format_dsc {
    container "${params.container__pandas}"
    publishDir "${params.output}/7_filtered_SSC/${specimen}/", mode: 'copy', overwrite: true
    
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
    publishDir "${params.output}/7_filtered_SSC/${specimen}/", mode: 'copy', overwrite: true
    
    input:
    tuple val(specimen), path("DSC.bam")
    path ref

    output:
    file "DSC.vcf.gz"

    script:
    template 'format_vcf.sh'

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

    // emit:
    // adduct_summary = call_adducts

}