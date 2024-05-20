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
    # We are at /docker/cellranger_pipeline
fi

# We are at /example, but works in any directory
echo "+----------------------------------------+"
echo "| Setting and creating workdir           |"
echo "+----------------------------------------+"
echo ""
WORKDIR=$(pwd)/example-data
mkdir -p $WORKDIR

echo "+----------------------------------------+"
echo "| Donwloading count test data            |"
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

if [[ "$1" == "local" ]]; then
    echo "+----------------------------------------+"
    echo "| Downloading count reference data       |"
    echo "+----------------------------------------+"
    echo ""
    if [[ ! -f $WORKDIR/refdata-gex-GRCh38-2020-A.tar.gz ]]; then
        wget -q https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz -O $WORKDIR/refdata-gex-GRCh38-2020-A.tar.gz
    fi

    echo "+----------------------------------------+"
    echo "| Extracting reference data              |"
    echo "+----------------------------------------+"
    tar -zxvf $WORKDIR/refdata-gex-GRCh38-2020-A.tar.gz -C $WORKDIR/

    echo "+----------------------------------------+"
    echo "| Downloading vdj test data              |"
    echo "+----------------------------------------+"
    echo ""
    if [[ ! -f $WORKDIR/sc5p_v2_hs_B_1k_multi_5gex_b_Multiplex_fastqs.tar ]]; then
        wget -q https://cf.10xgenomics.com/samples/cell-vdj/6.0.0/sc5p_v2_hs_B_1k_multi_5gex_b_Multiplex/sc5p_v2_hs_B_1k_multi_5gex_b_Multiplex_fastqs.tar -O $WORKDIR/sc5p_v2_hs_B_1k_multi_5gex_b_Multiplex_fastqs.tar
    fi

    echo "+----------------------------------------+"
    echo "| Extracting vdj test data               |"
    echo "+----------------------------------------+"
    tar -xvf $WORKDIR/sc5p_v2_hs_B_1k_multi_5gex_b_Multiplex_fastqs.tar -C $WORKDIR/

    echo "+----------------------------------------+"
    echo "| Downloading vdj reference data         |"
    echo "+----------------------------------------+"
    echo ""
    if [[ ! -f $WORKDIR/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0.tar.gz ]]; then
        wget -q https://cf.10xgenomics.com/supp/cell-vdj/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0.tar.gz -O $WORKDIR/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0.tar.gz
    fi

    echo "+----------------------------------------+"
    echo "| Extracting vdj reference data          |"
    echo "+----------------------------------------+"
    tar -zxvf $WORKDIR/refdata-cellranger-vdj-GRCh38-alts-ensembl-5.0.0.tar.gz -C $WORKDIR/
fi
