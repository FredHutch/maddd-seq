#!/bin/bash

# User must provide the S3 prefix containing a set of FASTQ files
DATASET=$1

# By default, only N sequences will be downloaded
N=${MAX_SEQUENCES:-1000000}
let "L=$N*4"

echo "Getting data from $DATASET"

# Get the list of FASTQ files in this prefix
FASTQ_PREFIX_LIST=$(aws s3 ls $DATASET --recursive --human-readable | grep fastq.gz | sed 's/.* //')
N_FASTQ=$( echo $FASTQ_PREFIX_LIST | tr ' ' '\n' | grep -c fastq.gz )

echo "Found ${N_FASTQ} files to download"

# Figure out the bucket name
BUCKET=$(echo $DATASET | sed 's/s3:\/\///' | sed 's/\/.*//')

# Make a manifest
i=0
while [[ -s pubweb-manifest-$i.csv ]]; do let "i=$i+1"; done
MANIFEST="pubweb-manifest-$i.csv"
echo "Writing out manifest to $MANIFEST"

# Write the header
echo "specimen,R1,R2" > $MANIFEST


download_fastq(){
    PREFIX=$1
    LOCAL_FP=$2
    echo "Downloading $N reads from s3://$BUCKET/$PREFIX to $LOCAL_FP"
    aws s3 cp s3://$BUCKET/$PREFIX - | gunzip -c | head -n $L | gzip -c > $LOCAL_FP
}

# Iterate over the files
for R1 in $FASTQ_PREFIX_LIST; do

    # Only consider the _R1_ files
    if [[ "${R1}" == *"_R1_"* ]]; then

        # Figure out the R2
        R2="${R1/_R1_/_R2_}"

        # Figure out the local file paths
        FP1="fastq/$(echo $R1 | sed 's/.*\///')"
        FP2="fastq/$(echo $R2 | sed 's/.*\///')"

        # Download the files
        download_fastq $R2 $FP2
        download_fastq $R1 $FP1
        
        # Figure out the specimen
        SPECIMEN="$(echo $R1 | sed 's/_R1_.*//' | sed 's/.*\///')"

        # Write out to the manifest
        echo "$SPECIMEN,$FP1,$FP2" >> $MANIFEST

    fi

done