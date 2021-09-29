#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

workflow variants_wf{

    take:
    bam_ch

    main:
    
    // Call adducts from the SSC data
    call_adducts(ssc_ch)

    // Call variants from the DSC data
    call_variants(dsc_ch)

    emit:
    adduct_summary = call_adducts

}