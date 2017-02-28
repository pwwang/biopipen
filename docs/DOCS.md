<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [Documentation for bioprocs v0.0.1](#documentation-for-bioprocs-v001)
  - [WXS](#wxs)
  - [WEB](#web)
    - [pDownloadPost](#pdownloadpost)
      - [description](#description)
      - [input](#input)
      - [output](#output)
      - [args](#args)
      - [requires](#requires)
    - [pDownloadGet](#pdownloadget)
      - [description](#description-1)
      - [input](#input-1)
      - [args](#args-1)
      - [output](#output-1)
  - [TCGA](#tcga)
    - [pSample2SubmitterID](#psample2submitterid)
      - [description](#description-2)
      - [input](#input-2)
      - [output](#output-2)
    - [pConvertExpFiles2Matrix](#pconvertexpfiles2matrix)
      - [description](#description-3)
      - [input](#input-3)
      - [output](#output-3)
      - [requires](#requires-1)
    - [pConvertMutFiles2Matrix](#pconvertmutfiles2matrix)
      - [description](#description-4)
      - [input](#input-4)
      - [output](#output-4)
  - [CHIPSEQ](#chipseq)
    - [pPeakToRegPotential](#ppeaktoregpotential)
      - [description](#description-5)
      - [input](#input-5)
      - [output](#output-5)
      - [args](#args-2)
  - [GSEA](#gsea)
    - [pMTarget2GTargetMat](#pmtarget2gtargetmat)
      - [description](#description-6)
      - [input](#input-6)
      - [output](#output-6)
      - [args](#args-3)
      - [requires](#requires-2)
    - [pSSGSEA](#pssgsea)
      - [description](#description-7)
      - [input](#input-7)
      - [output](#output-7)
      - [args](#args-4)
      - [requires](#requires-3)
  - [DEG](#deg)
    - [pCallByLimmaFromMatrix](#pcallbylimmafrommatrix)
      - [description](#description-8)
      - [input](#input-8)
      - [output](#output-8)
      - [args](#args-5)
      - [requires](#requires-4)
    - [pCallByLimmaFromFiles](#pcallbylimmafromfiles)
      - [description](#description-9)
      - [input](#input-9)
      - [output](#output-9)
      - [args](#args-6)
      - [requires](#requires-5)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->


# Documentation for bioprocs v0.0.1
A set of procs for bioinformatics using [pyppl](https://github.com/pwwang/pyppl)

## WXS

## WEB

###  pDownloadPost
#### description
- Download results by submitting a form, supporting pagination.

#### input
- `url`: the URL contains the form
- `submitbtn`: the submit button to click to submit the form (use Xpath).
- `nextbtn`: the button for next page (use Xpath)
- `params`: the params used to fill the form (JSON string or transformed from dict by json.dumps).

#### output
- `outdir:file`: The directory saves the results

#### args
- `interval`: seconds to wait between fetching each page. Default: 1

#### requires
- - [`Splinter`](https://splinter.readthedocs.io/en/latest/index.html)
- - [`Phantomjs`](http://phantomjs.org/)


###  pDownloadGet
#### description
- Download results by urls.

#### input
- `url`: the URLs to download

#### args
- `keepname`: bool, whether to keep the basename, otherwise use {{#}}.<ext>, default: True

#### output
- `outdir:file`: The directory saves the results


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


## CHIPSEQ

###  pPeakToRegPotential
#### description
- Convert peaks to regulatory potential score for each gene
- The formula is:
- ```
- -(0.5 + 4*di/d0)
- PC = sum (pi * e                  )
- ```
- Ref: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4489297/

#### input
- `peakfile:file`: The BED/peak file for peaks
- `genefile:file`: The BED file for gene coordinates

#### output
- `outfile:file`: The regulatory potential file for each gene

#### args
- `intensity`: `pi` in the formula. Boolean value, whether use the peak intensity or not, default: `True`,
- `geneformat`: The format for `genefile`, default: `ucsc+gz`. It could be:
- - ucsc or ucsc+gz: typically, you can download from http://hgdownload.cse.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz
- - bed or bed+gz: [format](https://genome.ucsc.edu/FAQ/FAQformat#format1), 4th column required as gene identity.
- `peakformat`: The format for `peakfile`, default: `peak`. It could be:
- - peak or peak+gz: (either [narrowPeak](https://genome.ucsc.edu/FAQ/FAQformat.html#format12) or [broadPeak](https://genome.ucsc.edu/FAQ/FAQformat.html#format13), the 7th column will be used as intensity
- - bed or bed+gz: [format](https://genome.ucsc.edu/FAQ/FAQformat#format1), 5th column will be used as intensity.
- `window`: `2 * d0` in the formula. The window where the peaks fall in will be consided, default: `100000`. 
- ```
- |--------- window ----------|
- |---- d0 -----|
- |--- 50K --- TSS --- 50K ---|
- ^ (peak center)
- |-- di --|
- ```


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

#### args
- `species`: The species used to convert gene names, default: human

#### requires
- [python-mygene](https://pypi.python.org/pypi/mygene/3.0.0)


###  pSSGSEA
#### description
- Single sample GSEA
- Refer to http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GCT for file format

#### input
- `gctfile:file`: the expression file
- `gmtfile:file`: the gmtfile for gene sets

#### output
- `outdir:file`: the output directory
- - `report.txt`: the enrichment report for each Gene set.
- - `RES_<GeneSet>.png`: the running ES plot for <GeneSet>
- - `normP_<GeneSet>.png`: the norminal P value plot for <GeneSet>

#### args
- `weightexp`: Exponential weight employed in calculation of enrichment scores. Default: 0.75
- `padjust`:   P value adjustment method, default 'bonferroni'. Can be "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
- `nperm`:     Number of permutations. Default: 10000

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

