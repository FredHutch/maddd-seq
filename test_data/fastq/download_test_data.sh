#!/bin/bash

set -Eeuo pipefail

N_READS=${N_READS:-100000}
SRA_LIST="SRR7897095 SRR7897096 SRR7897097 SRR7897100 SRR7897098"

download_sra(){
    SRA=$1

    fastq-dump -X $N_READS --split-files $SRA 
    gzip $SRA*.fastq

}

for SRA in $SRA_LIST; do

    echo "Downloading $SRA"
    download_sra $SRA

done