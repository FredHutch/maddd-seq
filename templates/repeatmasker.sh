#!/bin/bash

set -e

echo "Running RepeatMasker on ${genome_fasta}"

RepeatMasker \
    -qq \
    -pa ${task.cpus} \
    -species "${params.repeatmasker_species}" \
    "${genome_fasta}"

echo DONE