#!/bin/bash

set -euo pipefail

samtools \
    merge \
    --write-index \
    --threads ${task.cpus} \
    ${specimen}.${params.file_label}.bam \
    input_bam/*.bam
