#!/usr/bin/env bash

# This script generates the example data for the cnvkit-pipeline pipeline.

set -e

# Add an argument to indicate whether this is running locally
# or in the container
if [[ "$1" == "local" ]]; then
    echo "+----------------------------------------+"
    echo "| Running in local mode                  |"
    echo "+----------------------------------------+"
    echo ""
    # We are at /docker/pipseeker_pipeline
fi

# We are at /example, but works in any directory
echo "+----------------------------------------+"
echo "| Setting and creating workdir           |"
echo "+----------------------------------------+"
echo ""
WORKDIR=$(pwd)/example-data
mkdir -p $WORKDIR

echo "+----------------------------------------+"
echo "| Donwloading test data                  |"
echo "+----------------------------------------+"
echo ""
for tfile in \
    Sample1_Sub_R1.fastq.gz \
    Sample1_Sub_R2.fastq.gz \
    Sample2_Sub_R1.fastq.gz \
    Sample2_Sub_R2.fastq.gz;
do
    if [[ ! -f $WORKDIR/$tfile ]]; then
        wget -q https://github.com/noamteyssier/pipseq_nextflow/raw/refs/heads/main/data/sequences/$tfile -O "$WORKDIR/$(basename $tfile)"
    fi
done

if [[ "$1" == "local" ]]; then
    echo "+----------------------------------------+"
    echo "| Downloading reference data             |"
    echo "+----------------------------------------+"
    echo ""
    if [[ ! -f $WORKDIR/refdata-gex-GRCh38-2020-A.tar.gz ]]; then
        wget -q https://github.com/latchbio-workflows/wf-fluentbio-pipseeker/raw/refs/heads/main/unit_tests/test_data/STAR_test_index.zip -O $WORKDIR/refdata.tar.gz
    fi

    echo "+----------------------------------------+"
    echo "| Extracting reference data              |"
    echo "+----------------------------------------+"
    unzip -q $WORKDIR/refdata.tar.gz -d $WORKDIR/
fi
