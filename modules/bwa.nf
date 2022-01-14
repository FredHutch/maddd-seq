// Align reads with BWA MEM
process bwa {
    container "${params.container__bwa}"
    label "cpu_medium"
    
    input:
    tuple val(specimen), val(shard_ix), path(R1), path(R2)
    path ref

    output:
    tuple val(specimen), val(shard_ix), path("aligned.bam"), emit: bam

    script:
    template 'bwa.sh'

}

// Merge a collection of alignments in sorted BAM format
process merge_bam {
    container "${params.container__bwa}"
    publishDir "${params.output}/${params.subfolder}/", mode: 'copy', overwrite: true, enabled: params.publish
    label "cpu_medium"
    
    input:
    tuple val(specimen), path("input_bam/*.bam")

    output:
    path "${specimen}.${params.file_label}.*"

    script:
    template 'merge_bam.sh'

}
