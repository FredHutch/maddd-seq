#!/bin/bash

set -Eeuo pipefail

echo "Indexing $genome_fasta"

bwa index "${genome_fasta}"

echo "DONE"