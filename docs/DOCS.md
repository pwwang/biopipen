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
  - [ALGORITHM](#algorithm)
    - [pRWR](#prwr)
      - [description](#description-5)
      - [input](#input-5)
      - [output](#output-5)
      - [args](#args-2)
      - [requires](#requires-2)
  - [COMMON](#common)
    - [pSort](#psort)
      - [description](#description-6)
      - [input](#input-6)
      - [output](#output-6)
      - [args](#args-3)
  - [CHIPSEQ](#chipseq)
    - [pPeakToRegPotential](#ppeaktoregpotential)
      - [description](#description-7)
      - [input](#input-7)
      - [output](#output-7)
      - [args](#args-4)
  - [GSEA](#gsea)
    - [pMTarget2GTargetMat](#pmtarget2gtargetmat)
      - [description](#description-8)
      - [input](#input-8)
      - [output](#output-8)
      - [args](#args-5)
      - [requires](#requires-3)
    - [pIntersectGMT](#pintersectgmt)
      - [description](#description-9)
      - [input](#input-9)
      - [output](#output-9)
      - [args](#args-6)
      - [requires](#requires-4)
    - [pUnionGMT](#puniongmt)
      - [description](#description-10)
      - [input](#input-10)
      - [output](#output-10)
      - [args](#args-7)
      - [requires](#requires-5)
    - [pSSGSEA](#pssgsea)
      - [description](#description-11)
      - [input](#input-11)
      - [output](#output-11)
      - [args](#args-8)
      - [requires](#requires-6)
  - [SNPARRAY](#snparray)
    - [pSNP6Genotype](#psnp6genotype)
      - [description](#description-12)
      - [input](#input-12)
      - [output](#output-12)
      - [requires](#requires-7)
    - [pGenoToAvInput](#pgenotoavinput)
      - [description](#description-13)
      - [input](#input-13)
      - [output](#output-13)
      - [requires](#requires-8)
  - [DEG](#deg)
    - [pCallByLimmaFromMatrix](#pcallbylimmafrommatrix)
      - [description](#description-14)
      - [input](#input-14)
      - [output](#output-14)
      - [args](#args-9)
      - [requires](#requires-9)
    - [pCallByLimmaFromFiles](#pcallbylimmafromfiles)
      - [description](#description-15)
      - [input](#input-15)
      - [output](#output-15)
      - [args](#args-10)
      - [requires](#requires-10)

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
- `outfile:file`: The output file


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


## ALGORITHM

###  pRWR
#### description
- Do random walk with restart (RWR)

#### input
- `Wfile:file`: The adjecent matrix
- `Efile:file`: The start vector

#### output
- `outfile:file`: The output of final probabilities

#### args
- `c`:       The restart probability
- `eps`:     The convergent cutoff || R(i+1) - R(i) ||
- `tmax`:    Max iterations to stop
- `Wformat`: The format of Wfile, rds or mat/txt, default: rds
- `Eformat`: The format of Efile, rds or mat/txt, default: rds
- `Rformat`: The format of the output file, rds or mat/txt, default: rds
- `normW`:   Weather to normalize W or not, default False. Laplacian normalization is used (more to add).
- `normE`:   Weather to normalize E or not, default False. E will be normalized as: E = E/sum(E)

#### requires
- if normW = True, R package `NetPreProc` is required.


## COMMON

###  pSort
#### description
- Sort file using linux command `sort`

#### input
- `infile:file`: The input file

#### output
- `outfile:file`: The output file

#### args
- `sort-args`: The arguments used by `sort`


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


###  pIntersectGMT
#### description
- Get the intersect gene set from multiple gmt files
- To do intersect for more than 2 files: gmtfile1, gmtfile2, gmtfile3:
- ```
- pIntersectGMT.input = {pIntersectGMT.input: channel.create([(gmtfile1, gmtfile2)])}
- pIntersectGMT2 = pIntersectGMT.copy()
- pIntersectGMT2.depends = pIntersectGMT
- pIntersectGMT2.input   = {pIntersectGMT2.input.keys()[0]: lambda ch: ch.insert(0, gmtfile3)}
- ```

#### input
- `gmtfile1:file`: The 1st gmt file
- `gmtfile2:file`: The 2nd gmt file

#### output
- `outdir:file`: the output gmtfile

#### args
- `geneformat`: The gene names in gene set. Default: "symbol,alias". Available values see mygene docs.
- `gz`: whether the files are with gz format, default: False. If `gz = True`, output file will be also gzipped.
- `species`: The species, used for gene name conversion in mygene, default: human

#### requires
- [python-mygene](https://pypi.python.org/pypi/mygene/3.0.0) 


###  pUnionGMT
#### description
- Get the union gene set from multiple gmt files
- To do union for more than 2 files: gmtfile1, gmtfile2, gmtfile3:
- ```
- pUnionGMT.input = {pIntersectGMT.input: channel.create([(gmtfile1, gmtfile2)])}
- pUnionGMT2 = pIntersectGMT.copy()
- pUnionGMT2.depends = pIntersectGMT
- pUnionGMT2.input   = {pUnionGMT.input.keys()[0]: lambda ch: ch.insert(0, gmtfile3)}
- ```

#### input
- `gmtfile1:file`: The 1st gmt file
- `gmtfile2:file`: The 2nd gmt file

#### output
- `outdir:file`: the output gmtfile

#### args
- `geneformat`: The gene names in gene set. Default: "symbol,alias". Available values see mygene docs.
- `gz`: whether the files are with gz format, default: False. If `gz = True`, output file will be also gzipped.
- `species`: The species, used for gene name conversion in mygene, default: human

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


## SNPARRAY

###  pSNP6Genotype
#### description
- Call genotypes from GenomeWideSNP_6 CEL file

#### input
- `celfile:file`: the CEL file

#### output
- `outfile:file`: the outfile containing probe name and genotypes
- format: <Probe name>\t<genotype>
- <genotype> = 0: AA, 1: AB, 2: BB

#### requires
- [bioconductor-crlmm](http://bioconductor.org/packages/release/bioc/html/crlmm.html)


###  pGenoToAvInput
#### description
- Convert the genotype called by pSNP6Genotype to [ANNOVAR input file](http://annovar.openbioinformatics.org/en/latest/user-guide/input/#annovar-input-file) using dbSNP identifiers.	

#### input
- `genofile:file`: the genofile generated by pSNP6Genotype, must be sorted by probe names
- `annofile:flie`: the annotation file downloaded from http://www.affymetrix.com/support/technical/annotationfilesmain.affx
- Could be in .gz format

#### output
- `outfile:file`: the avinput file

#### requires
- [python-read2](https://github.com/pwwang/read2)


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
- `paired`: whether the samples are paired

#### requires
- [limma](https://bioconductor.org/packages/release/bioc/html/limma.html)

