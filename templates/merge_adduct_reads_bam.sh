#!/bin/bash

set -euo pipefail

samtools \
    merge \
    --write-index \
    --threads ${task.cpus} \
    ${filtering}.adduct_reads.bam \
    input_bam/*.bam

samtools index ${filtering}.adduct_reads.bam