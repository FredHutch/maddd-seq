#!/usr/bin/env nextflow

// Using DSL-2
nextflow.enable.dsl=2

workflow manifest_wf{

    take:
    manifest_fp

    main:

    reads_ch = Channel.from(
        file(manifest_fp)
    ).splitCsv(
        header: true
    ).map {
        r -> [
            r["specimen"], 
            file(r["R1"]),
            file(r["R2"])
        ]
    }

    emit:
    reads = reads_ch

}