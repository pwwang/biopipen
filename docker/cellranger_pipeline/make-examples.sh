#!/usr/bin/env bash

# This script generates the example data for the cnvkit-pipeline pipeline.

set -e

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
    Sample_X_S1_L001_R1_001.fastq.gz \
    Sample_X_S1_L001_R2_001.fastq.gz \
    Sample_Y_S1_L001_R1_001.fastq.gz \
    Sample_Y_S1_L001_R2_001.fastq.gz \
    Sample_Y_S1_L002_R1_001.fastq.gz \
    Sample_Y_S1_L002_R2_001.fastq.gz;
do
    if [[ ! -f $WORKDIR/$tfile ]]; then
        wget -q https://raw.githubusercontent.com/nf-core/test-datasets/scrnaseq/testdata/cellranger/$tfile -O "$WORKDIR/$(basename $tfile)"
    fi
done

# echo "+----------------------------------------+"
# echo "| Downloading reference data             |"
# echo "+----------------------------------------+"
# echo ""
# if [[ ! -f $WORKDIR/refdata-gex-GRCh38-2020-A.tar.gz ]]; then
#     wget -q wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz -O $WORKDIR/refdata-gex-GRCh38-2020-A.tar.gz
# fi

# echo "+----------------------------------------+"
# echo "| Extracting reference data              |"
# echo "+----------------------------------------+"
# tar -zxvf $WORKDIR/refdata-gex-GRCh38-2020-A.tar.gz -C $WORKDIR/
