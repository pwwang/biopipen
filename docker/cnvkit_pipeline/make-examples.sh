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
    genome/genome.fasta \
    illumina/bam/test.paired_end.sorted.bam \
    illumina/bam/test2.paired_end.sorted.bam;
do
    if [[ ! -f $WORKDIR/$tfile ]]; then
        wget -q https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/genomics/homo_sapiens/$tfile -O $WORKDIR/$(basename $tfile)
    fi
done

echo "+----------------------------------------+"
echo "| Donwloading refFlat                    |"
echo "+----------------------------------------+"
echo ""
if [[ ! -f $WORKDIR/refFlat.txt ]]; then
    wget -q https://hgdownload.cse.ucsc.edu/goldenpath/hg19/database/refFlat.txt.gz -O $WORKDIR/refFlat.txt.gz
    gunzip $WORKDIR/refFlat.txt.gz
fi

echo "+----------------------------------------+"
echo "| Making the metadata file               |"
echo "+----------------------------------------+"
echo ""

metafile=$WORKDIR/metadata.txt
echo -e "Sample\tBam\tGuessBaits" > $metafile
echo -e "Test\t$WORKDIR/test.paired_end.sorted.bam\t1" >> $metafile
echo -e "Test2\t$WORKDIR/test2.paired_end.sorted.bam\t1" >> $metafile

echo "+----------------------------------------+"
echo "| Preparing the data for pipen-board     |"
echo "+----------------------------------------+"
echo ""
# Only works in the docker container
if [ -d /biopipen ]; then
    # L3dvcmtkaXIvLnBpcGVu is the base64 encoded string of "/workdir/.pipen"
    cp /biopipen/docker/cnvkit_pipeline/example.json \
        /biopipen/.pipen-board/biopipen-ns-cnvkit-pipeline-cnvkitpipeline.Example.L3dvcmtkaXIvLnBpcGVu.json
else
    mkdir -p ~/.pipen-board
    wget https://raw.githubusercontent.com/pwwang/biopipen/master/docker/cnvkit_pipeline/example.json -O ~/.pipen-board/biopipen-ns-cnvkit-pipeline-cnvkitpipeline.Example.L3dvcmtkaXIvLnBpcGVu.json
fi
