#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

workflow manifest_wf{

    take:
    manifest_fp

    main:
    parse_manifest(manifest_fp)

    emit:
    reads = parse_manifest.out

}