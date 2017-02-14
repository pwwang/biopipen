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
  - [DEG](#deg)
    - [pCallByLimmaFromMatrix](#pcallbylimmafrommatrix)
      - [description](#description-2)
      - [input](#input-2)
      - [output](#output-2)
      - [args](#args)
      - [requires](#requires-1)
    - [pCallByLimmaFromFiles](#pcallbylimmafromfiles)
      - [description](#description-3)
      - [input](#input-3)
      - [output](#output-3)
      - [args](#args-1)
      - [requires](#requires-2)

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

