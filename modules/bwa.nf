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