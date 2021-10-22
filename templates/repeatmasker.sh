#!/bin/bash

set -e

echo "Running RepeatMasker on ${genome_fasta}"

RepeatMasker -species "${params.repeatmasker_species}" "${genome_fasta}"

echo DONE