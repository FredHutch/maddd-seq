#!/bin/bash

set -euo pipefail

samtools merge -o ${specimen}.${params.file_label}.bam input_bam/*.bam

samtools index ${specimen}.${params.file_label}.bam
