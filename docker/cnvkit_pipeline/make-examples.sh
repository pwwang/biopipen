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
    reference/human_g1k_v37_decoy.small.fasta \
    testdata/target.bed \
    testdata/recalbam/9876T.recal.bam \
    testdata/recalbam/9876T.recal.bai \
    testdata/recalbam/1234N.recal.bam \
    testdata/recalbam/1234N.recal.bai;
do
    if [[ ! -f $WORKDIR/$tfile ]]; then
        wget -q https://raw.githubusercontent.com/nf-core/test-datasets/sarek/$tfile -O $WORKDIR/$(basename $tfile)
    fi
done

echo "+----------------------------------------+"
echo "| Donwloading refFlat                    |"
echo "+----------------------------------------+"
echo ""
if [[ ! -f $WORKDIR/refFlat.txt ]]; then
    wget -q https://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/refFlat.txt.gz -O $WORKDIR/refFloat.txt.gz
    gunzip $WORKDIR/refFlat.txt.gz
fi

echo "+----------------------------------------+"
echo "| Making the metadata file               |"
echo "+----------------------------------------+"
echo ""

metafile=$WORKDIR/metadata.txt
echo -e "Sample\tBam\tGroup" > $metafile
echo -e "9876T\t$WORKDIR/9876T.recal.bam\tTumor" >> $metafile
echo -e "1234N\t$WORKDIR/1234N.recal.bam\tNormal" >> $metafile

echo "+----------------------------------------+"
echo "| Preparing the data for pipen-board     |"
echo "+----------------------------------------+"
echo ""
# Only works in the docker container
if [ -d /biopipen ]; then
    cp /biopipen/docker/cnvkit_pipeline/example.json \
        /biopipen/.pipen-board/biopipen-ns-cnvkit-pipeline-cnvkitpipeline.Example.0000-00-00_00-00-00.json
else
    mkdir -p ~/.pipen-board
    wget https://raw.githubusercontent.com/pwwang/biopipen/master/docker/cnvkit_pipeline/example.json -O ~/.pipen-board/biopipen-ns-cnvkit-pipeline-cnvkitpipeline.Example.0000-00-00_00-00-00.json
fi
