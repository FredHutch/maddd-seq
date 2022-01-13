#!/bin/bash

set -euo pipefail

cat R1/* > ${specimen}_R1.fastq.gz
cat R2/* > ${specimen}_R2.fastq.gz
