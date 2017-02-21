<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [Documentation for bioprocs v0.0.1](#documentation-for-bioprocs-v001)
  - [TCGA](#tcga)
    - [pSample2SubmitterID](#psample2submitterid)
      - [description](#description)
      - [input](#input)
      - [output](#output)
    - [pConvertExpFiles2Matrix](#pconvertexpfiles2matrix)
      - [description](#description-1)
      - [input](#input-1)
      - [output](#output-1)
      - [requires](#requires)
    - [pConvertMutFiles2Matrix](#pconvertmutfiles2matrix)
      - [description](#description-2)
      - [input](#input-2)
      - [output](#output-2)
  - [GSEA](#gsea)
    - [pMTarget2GTargetMat](#pmtarget2gtargetmat)
      - [description](#description-3)
      - [input](#input-3)
      - [output](#output-3)
      - [requires](#requires-1)
  - [DEG](#deg)
    - [pCallByLimmaFromMatrix](#pcallbylimmafrommatrix)
      - [description](#description-4)
      - [input](#input-4)
      - [output](#output-4)
      - [args](#args)
      - [requires](#requires-2)
    - [pCallByLimmaFromFiles](#pcallbylimmafromfiles)
      - [description](#description-5)
      - [input](#input-5)
      - [output](#output-5)
      - [args](#args-1)
      - [requires](#requires-3)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->


# Documentation for bioprocs v0.0.1
A set of procs for bioinformatics using [pyppl](https://github.com/pwwang/pyppl)

## TCGA

###  pSample2SubmitterID
#### description
- convert TCGA sample names with submitter id with metadata and sample containing folder

#### input
- `dir:file`:    the directory containing the samples
- `mdfile:file`: the metadata file

#### output
- `outdir:file`: the directory containing submitter-id named files


###  pConvertExpFiles2Matrix
#### description
- convert TCGA expression files to expression matrix, and convert sample name to submitter id

#### input
- `dir:file`:    the directory containing the samples
- `mdfile:file`: the metadata file

#### output
- `outfile:file`:the output matrix

#### requires
- [python-mygene](https://pypi.python.org/pypi/mygene/3.0.0)


###  pConvertMutFiles2Matrix
#### description
- convert TCGA mutation files (vcf.gz) to mut matrix, and convert sample name to submitter id

#### input
- `dir:file`:    the directory containing the samples
- `mdfile:file`: the metadata file

#### output
- `outfile:file`:the output matrix


## GSEA

###  pMTarget2GTargetMat
#### description
- Convert motif target from MSigDB database (i.e. c3.tft.v5.2.entrez.gmt from GSEA to gene-target matrix
- You also have to have a map of motif name to genes (https://github.com/andrewdyates/transcription_factors/blob/master/gsea_msigdb/transfac_id_to_genes_raw.tab)

#### input
- `gmtfile:file`: typically c3.tft.v5.2.entrez.gmt (have to be in entrez format)
- `mapfile:file`: the motif-gene name map file

#### output
- `outfile:file`: the gene-target matrix

#### requires
- [python-mygene](https://pypi.python.org/pypi/mygene/3.0.0)


## DEG

###  pCallByLimmaFromMatrix
#### description
- Call DEG from expressoin matrix, where column names must in accordant order of <group>

#### input
- `matfile:file`: the expression matrix
- `group1`:       columns of group1 (separated by comma)
- `group2`:       columns of group2 (separated by comma)
- `group1name`:   the name of group1
- `group2name`:   the name of group2

#### output
- degfile:file: the output file containing DEGs

#### args
- `pval`: the cutoff of DEGs (default: .05)

#### requires
- [limma](https://bioconductor.org/packages/release/bioc/html/limma.html)


###  pCallByLimmaFromFiles
#### description
- Call DEG from expression files

#### input
- `expdir:file`:  the directory containing expression files
- `group1`:       columns of group1 (separated by comma)
- `group2`:       columns of group2 (separated by comma)
- `group1name`:   the name of group1
- `group2name`:   the name of group2   

#### output
- `degfile:file`: the output file containing DEGs

#### args
- `pval`: the cutoff of DEGs (default: .05)

#### requires
- [limma](https://bioconductor.org/packages/release/bioc/html/limma.html)

