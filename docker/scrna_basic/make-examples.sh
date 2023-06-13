#!/usr/bin/env bash

# This script generates the example data for the scrna-basic pipeline.
# It is intended to be run inside the docker container.

set -e

# We are at /example, but works in any directory
echo "+----------------------------------------+"
echo "| Setting workdir to current directory   |"
echo "+----------------------------------------+"
echo ""
WORKDIR=$(pwd)

echo "+----------------------------------------+"
echo "| Donwloading example data               |"
echo "+----------------------------------------+"
echo ""
if [[ ! -f $WORKDIR/ifnb.rda ]]; then
    wget https://seurat.nygenome.org/src/contrib/ifnb.SeuratData_3.1.0.tar.gz -O $WORKDIR/ifnb.SeuratData_3.1.0.tar.gz
    tar zxvf ifnb.SeuratData_3.1.0.tar.gz -C $WORKDIR --strip-components=2 ifnb.SeuratData/data/ifnb.rda
    rm -f $WORKDIR/ifnb.SeuratData_3.1.0.tar.gz
fi

echo "+----------------------------------------+"
echo "| Downloading ref data for clustering    |"
echo "+----------------------------------------+"
echo ""
if [[ ! -f $WORKDIR/pbmc_multimodal.h5seurat ]]; then
    wget https://atlas.fredhutch.org/data/nygc/multimodal/pbmc_multimodal.h5seurat -O $WORKDIR/pbmc_multimodal.h5seurat
fi

echo "+----------------------------------------+"
echo "| Downloading KEGG pathways              |"
echo "+----------------------------------------+"
echo ""
if [[ ! -f $WORKDIR/KEGG_pathways.gmt ]]; then
    wget https://data.broadinstitute.org/gsea-msigdb/msigdb/release/7.5.1/c2.cp.kegg.v7.5.1.symbols.gmt -O $WORKDIR/KEGG_pathways.gmt
fi

echo "+----------------------------------------+"
echo "| Downloading sctype db                  |"
echo "+----------------------------------------+"
echo ""
if [[ ! -f $WORKDIR/ScTypeDB_full.xlsx ]]; then
    wget https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx -O $WORKDIR/ScTypeDB_full.xlsx
fi

echo "+----------------------------------------+"
echo "| Creating data directories              |"
echo "+----------------------------------------+"
echo ""
mkdir -p $WORKDIR/data/CTRL  # /example/data/CTRL
mkdir -p $WORKDIR/data/STIM # /example/data/STIM

echo "+----------------------------------------+"
echo "| Formatting the data                    |"
echo "+----------------------------------------+"
echo ""

echo "
library(Seurat)
library(DropletUtils)

set.seed(123)
setwd('$WORKDIR')

# Load the data
load('ifnb.rda')

# Split the data into two
objs <- SplitObject(ifnb, split.by = 'stim')

counts_ctrl <- GetAssayData(objs\$CTRL, slot = 'counts')
counts_stim <- GetAssayData(objs\$STIM, slot = 'counts')

# Save the data
write10xCounts('data/CTRL', counts_ctrl, version = '3', overwrite = TRUE)
write10xCounts('data/STIM', counts_stim, version = '3', overwrite = TRUE)
" | R --no-save

echo "+----------------------------------------+"
echo "| Generating the sample file             |"
echo "+----------------------------------------+"
echo ""
sample_file="$WORKDIR/example.txt"
echo -e "Sample\tRNADir" > $sample_file
echo -e "CTRL\t$WORKDIR/data/CTRL" >> $sample_file
echo -e "STIM\t$WORKDIR/data/STIM" >> $sample_file

echo "+----------------------------------------+"
echo "| Preparing the data for pipen-board     |"
echo "+----------------------------------------+"
echo ""
mkdir -p ~/.pipen-board/
cp /biopipen/docker/example.json biopipen-ns-scrna-basic-scrnabasic.Example.1.json
