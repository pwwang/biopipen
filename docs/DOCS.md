<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [Documentation for bioprocs v0.0.1](#documentation-for-bioprocs-v001)
  - [WXS](#wxs)
    - [pTrimmomaticPE](#ptrimmomaticpe)
      - [description](#description)
      - [input](#input)
      - [output](#output)
      - [args](#args)
      - [requires](#requires)
    - [pTrimmomaticSE](#ptrimmomaticse)
      - [description](#description-1)
      - [input](#input-1)
      - [output](#output-1)
      - [args](#args-1)
      - [requires](#requires-1)
    - [pAlignPEByBWA](#palignpebybwa)
      - [description](#description-2)
      - [input](#input-2)
      - [output](#output-2)
      - [args](#args-2)
      - [requires](#requires-2)
    - [pAlignSEByBWA](#palignsebybwa)
      - [description](#description-3)
      - [input](#input-3)
      - [output](#output-3)
      - [args](#args-3)
      - [requires](#requires-3)
    - [pSortSam](#psortsam)
      - [description](#description-4)
      - [input](#input-4)
      - [output](#output-4)
      - [args](#args-4)
      - [requires](#requires-4)
    - [pMarkDup](#pmarkdup)
      - [description](#description-5)
      - [input](#input-5)
      - [output](#output-5)
      - [args](#args-5)
      - [requires](#requires-5)
    - [pIndexBam](#pindexbam)
      - [description](#description-6)
      - [input](#input-6)
      - [output](#output-6)
      - [args](#args-6)
      - [requires](#requires-6)
    - [pCNVnator](#pcnvnator)
      - [description](#description-7)
      - [input](#input-7)
      - [output](#output-7)
      - [args](#args-7)
      - [requires](#requires-7)
  - [WEB](#web)
    - [pDownloadPost](#pdownloadpost)
      - [description](#description-8)
      - [input](#input-8)
      - [output](#output-8)
      - [args](#args-8)
      - [requires](#requires-8)
    - [pDownloadGet](#pdownloadget)
      - [description](#description-9)
      - [input](#input-9)
      - [args](#args-9)
      - [output](#output-9)
  - [TCGA](#tcga)
    - [pSample2SubmitterID](#psample2submitterid)
      - [description](#description-10)
      - [input](#input-10)
      - [output](#output-10)
    - [pConvertExpFiles2Matrix](#pconvertexpfiles2matrix)
      - [description](#description-11)
      - [input](#input-11)
      - [output](#output-11)
      - [requires](#requires-9)
    - [pConvertMutFiles2Matrix](#pconvertmutfiles2matrix)
      - [description](#description-12)
      - [input](#input-12)
      - [output](#output-12)
  - [ALGORITHM](#algorithm)
    - [pRWR](#prwr)
      - [description](#description-13)
      - [input](#input-13)
      - [output](#output-13)
      - [args](#args-10)
      - [requires](#requires-10)
  - [GATK](#gatk)
    - [pRealignerTargetCreator](#prealignertargetcreator)
      - [description](#description-14)
      - [input](#input-14)
      - [output](#output-14)
      - [args](#args-11)
      - [requires](#requires-11)
    - [pIndelRealigner](#pindelrealigner)
      - [description](#description-15)
      - [input](#input-15)
      - [output](#output-15)
      - [args](#args-12)
      - [requires](#requires-12)
    - [pBaseRecalibrator](#pbaserecalibrator)
      - [description](#description-16)
      - [input](#input-16)
      - [output](#output-16)
      - [args](#args-13)
      - [requires](#requires-13)
    - [pPrintReads](#pprintreads)
      - [description](#description-17)
      - [input](#input-17)
      - [output](#output-17)
      - [args](#args-14)
      - [requires](#requires-14)
    - [pHaplotypeCaller](#phaplotypecaller)
      - [description](#description-18)
      - [input](#input-18)
      - [output](#output-18)
      - [args](#args-15)
      - [requires](#requires-15)
    - [pSelectVariants](#pselectvariants)
      - [description](#description-19)
      - [input](#input-19)
      - [output](#output-19)
      - [args](#args-16)
      - [requires](#requires-16)
    - [pVariantFiltration](#pvariantfiltration)
      - [description](#description-20)
      - [input](#input-20)
      - [output](#output-20)
      - [args](#args-17)
      - [requires](#requires-17)
  - [COMMON](#common)
    - [pSort](#psort)
      - [description](#description-21)
      - [input](#input-21)
      - [output](#output-21)
      - [args](#args-18)
  - [CHIPSEQ](#chipseq)
    - [pPeakToRegPotential](#ppeaktoregpotential)
      - [description](#description-22)
      - [input](#input-22)
      - [output](#output-22)
      - [args](#args-19)
  - [GSEA](#gsea)
    - [pMTarget2GTargetMat](#pmtarget2gtargetmat)
      - [description](#description-23)
      - [input](#input-23)
      - [output](#output-23)
      - [args](#args-20)
      - [requires](#requires-18)
    - [pIntersectGMT](#pintersectgmt)
      - [description](#description-24)
      - [input](#input-24)
      - [output](#output-24)
      - [args](#args-21)
      - [requires](#requires-19)
    - [pUnionGMT](#puniongmt)
      - [description](#description-25)
      - [input](#input-25)
      - [output](#output-25)
      - [args](#args-22)
      - [requires](#requires-20)
    - [pSSGSEA](#pssgsea)
      - [description](#description-26)
      - [input](#input-26)
      - [output](#output-26)
      - [args](#args-23)
      - [requires](#requires-21)
  - [SNPARRAY](#snparray)
    - [pSNP6Genotype](#psnp6genotype)
      - [description](#description-27)
      - [input](#input-27)
      - [output](#output-27)
      - [requires](#requires-22)
    - [pGenoToAvInput](#pgenotoavinput)
      - [description](#description-28)
      - [input](#input-28)
      - [output](#output-28)
      - [requires](#requires-23)
  - [DEG](#deg)
    - [pCallByLimmaFromMatrix](#pcallbylimmafrommatrix)
      - [description](#description-29)
      - [input](#input-29)
      - [output](#output-29)
      - [args](#args-24)
      - [requires](#requires-24)
    - [pCallByLimmaFromFiles](#pcallbylimmafromfiles)
      - [description](#description-30)
      - [input](#input-30)
      - [output](#output-30)
      - [args](#args-25)
      - [requires](#requires-25)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->


# Documentation for bioprocs v0.0.1
A set of procs for bioinformatics using [pyppl](https://github.com/pwwang/pyppl)

## WXS

###  pTrimmomaticPE
#### description
- Trimming Illumina NGS paired-end data

#### input
- `fqfile1:file`: The 1st fastq file (could be in .gz format)
- `fqfile2:file`: The 2nd fastq file

#### output
- `outfile1:file`: The 1st output file
- `outfile2:file`: The 2nd output file

#### args
- `bin`:    The trimmomatic executable, default: "trimmomatic"
- `phred`:  "phred33" (default) or "phred64"
- `params`: Other params for trimmomatric, default: "ILLUMINACLIP:{adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
- - have to replace `{adapter}` with the path of the adapter file
- `nthread`: 1

#### requires
- [trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic)


###  pTrimmomaticSE
#### description
- Trimming Illumina NGS single-end data

#### input
- `fqfile:file`: The fastq file (could be in .gz format)

#### output
- `outfile:file`: The output file

#### args
- `bin`:    The trimmomatic executable, default: "trimmomatic"
- `phred`:  "phred33" (default) or "phred64"
- `params`: Other params for trimmomatric, default: "ILLUMINACLIP:{adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"
- - have to replace `{adapter}` with the path of the adapter file
- `nthread`: 1

#### requires
- [trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic)


###  pAlignPEByBWA
#### description
- Align paired-end reads to reference genome using bwa mem

#### input
- `infile1:file`: read file 1 (fastq, or fastq gzipped)
- `infile2:file`: read file 2 (fastq, or fastq gzipped)

#### output
- `outfile:file`: The output sam file

#### args
- `bin`:    The bwa executable, default: bwa
- `params`: Other params for bwa mem, default: "-M"
- `nthread`: 1
- `reffile`: The reference file

#### requires
- [bwa](https://github.com/lh3/bwa)


###  pAlignSEByBWA
#### description
- Align paired-end reads to reference genome using bwa mem

#### input
- `infile:file`:  read file (fastq, or fastq gzipped)

#### output
- `outfile:file`: The output sam file

#### args
- `bin`:    The bwa executable, default: bwa
- `params`: Other params for bwa mem, default: "-M"
- `nthread`: 1
- `reffile`: The reference file, required

#### requires
- [bwa](https://github.com/lh3/bwa)


###  pSortSam
#### description
- Use `picard SortSam` to sort sam or bam file

#### input
- `infile:file`:  The sam or bam file to be sorted

#### output
- `outfile:file`: The sorted sam or bam file

#### args
- `bin`:    The picard executable, default: "picard SortSam"
- `order`:  The sort order, default: coordinate. Possible: unsorted, queryname, coordinate, duplicate
- `outtype`:The type of output file, sam or bam. Default: bam

#### requires
- [picard](http://broadinstitute.github.io/picard/command-line-overview.html)


###  pMarkDup
#### description
- Use `picard MarkDuplicates` to  mark duplicates for bam file

#### input
- `infile:file`:  The bam file 

#### output
- `outfile:file`: The marked bam file

#### args
- `bin`:    The picard executable, default: "picard MarkDuplicates"
- `params`:  Other parameters for picard MarkDuplicates, default: ""

#### requires
- [picard](http://broadinstitute.github.io/picard/command-line-overview.html)


###  pIndexBam
#### description
- Use `picard BuildBamIndex` to index bam file

#### input
- `infile:file`:  The bam file 

#### output
- `outfile:file`: The same bam file (link) but with .bai file in `proc.outdir`

#### args
- `bin`:    The picard executable, default: "picard BuildBamIndex"
- `params`:  Other parameters for picard , default: ""

#### requires
- [picard](http://broadinstitute.github.io/picard/command-line-overview.html)


###  pCNVnator
#### description
- Use `CNVnator` to call CNVs from bam file

#### input
- `infile:file`:  The bam file 

#### output
- `outfile:file`: The vcf file

#### args
- `bin`:      The CNVnator executable, default: "cnvnator"
- `bin-vcf`:  The converter executable to convert CNVnator results to vcf, default: "cnvnator2VCF.pl"
- `binsize`:  The bin_size, default: 100
- `genome`:   The genome: default: hg19
- `chrom`:    Chromosome names, default: "" (all chromosomes)
- `chrdir`:   The dir contains reference sequence of chromosomes, default: "" (don't specify)
- 

#### requires
- [CNVnator](https://github.com/abyzovlab/CNVnator)


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
- [`Splinter`](https://splinter.readthedocs.io/en/latest/index.html)
- [`Phantomjs`](http://phantomjs.org/)


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


## GATK

###  pRealignerTargetCreator
#### description
- The local realignment process is designed to consume one or more BAM files and to locally realign reads such that the number of mismatching bases is minimized across all the reads. In general, a large percent of regions requiring local realignment are due to the presence of an insertion or deletion (indels) in the individual's genome with respect to the reference genome. Such alignment artifacts result in many bases mismatching the reference near the misalignment, which are easily mistaken as SNPs. Moreover, since read mapping algorithms operate on each read independently, it is impossible to place reads on the reference genome such that mismatches are minimized across all reads. Consequently, even when some reads are correctly mapped with indels, reads covering the indel near just the start or end of the read are often incorrectly mapped with respect the true indel, also requiring realignment. Local realignment serves to transform regions with misalignments due to indels into clean reads containing a consensus indel suitable for standard variant discovery approaches.
- Note that indel realignment is no longer necessary for variant discovery if you plan to use a variant caller that performs a haplotype assembly step, such as HaplotypeCaller or MuTect2. However it is still required when using legacy callers such as UnifiedGenotyper or the original MuTect. There are 2 steps to the realignment process:
- - Determining (small) suspicious intervals which are likely in need of realignment (RealignerTargetCreator)
- - Running the realigner over those intervals (see the IndelRealigner tool)
- For more details, see [the indel realignment method documentation](http://www.broadinstitute.org/gatk/guide/article?id=38).

#### input
- `infile:file`:  The aligned bam file 

#### output
- `outfile:file`: A list of target intervals to pass to the IndelRealigner.

#### args
- `bin`:     The gatk executable, default: "gatk -T RealignerTargetCreator"
- `params`:  Other parameters for RealignerTargetCreator, default: ""
- `reffile`: The reference file

#### requires
- [GATK](https://software.broadinstitute.org/gatk)


###  pIndelRealigner 
#### description
- The local realignment process is designed to consume one or more BAM files and to locally realign reads such that the number of mismatching bases is minimized across all the reads. In general, a large percent of regions requiring local realignment are due to the presence of an insertion or deletion (indels) in the individual's genome with respect to the reference genome. Such alignment artifacts result in many bases mismatching the reference near the misalignment, which are easily mistaken as SNPs. Moreover, since read mapping algorithms operate on each read independently, it is impossible to place reads on the reference genome such at mismatches are minimized across all reads. Consequently, even when some reads are correctly mapped with indels, reads covering the indel near just the start or end of the read are often incorrectly mapped with respect the true indel, also requiring realignment. Local realignment serves to transform regions with misalignments due to indels into clean reads containing a consensus indel suitable for standard variant discovery approaches.
- Note that indel realignment is no longer necessary for variant discovery if you plan to use a variant caller that performs a haplotype assembly step, such as HaplotypeCaller or MuTect2. However it is still required when using legacy callers such as UnifiedGenotyper or the original MuTect.
- There are 2 steps to the realignment process:
- - Determining (small) suspicious intervals which are likely in need of realignment (see the RealignerTargetCreator tool)
- - Running the realigner over those intervals (IndelRealigner)
- For more details, see [the indel realignment method documentation](http://www.broadinstitute.org/gatk/guide/article?id=38).

#### input
- `bamfile:file`: The aligned bam file
- `intfile:file`: Intervals file output from RealignerTargetCreator

#### output
- `outfile:file`: A realigned version of input BAM file.

#### args
- `bin`:     The gatk executable, default: "gatk -T IndelRealigner"
- `params`:  Other parameters for IndelRealigner, default: ""
- `reffile`: The reference file

#### requires
- [GATK](https://software.broadinstitute.org/gatk)


###  pBaseRecalibrator  
#### description
- Variant calling algorithms rely heavily on the quality scores assigned to the individual base calls in each sequence read. These scores are per-base estimates of error emitted by the sequencing machines. Unfortunately the scores produced by the machines are subject to various sources of systematic technical error, leading to over- or under-estimated base quality scores in the data. Base quality score recalibration (BQSR) is a process in which we apply machine learning to model these errors empirically and adjust the quality scores accordingly. This allows us to get more accurate base qualities, which in turn improves the accuracy of our variant calls. The base recalibration process involves two key steps: first the program builds a model of covariation based on the data and a set of known variants (which you can bootstrap if there is none available for your organism), then it adjusts the base quality scores in the data based on the model. There is an optional but highly recommended step that involves building a second model and generating before/after plots to visualize the effects of the recalibration process. This is useful for quality control purposes. This tool performs the first step described above: it builds the model of covariation and produces the recalibration table. It operates only at sites that are not in dbSNP; we assume that all reference mismatches we see are therefore errors and indicative of poor base quality. This tool generates tables based on various user-specified covariates (such as read group, reported quality score, cycle, and context). Assuming we are working with a large amount of data, we can then calculate an empirical probability of error given the particular covariates seen at this site, where p(error) = num mismatches / num observations. The output file is a table (of the several covariate values, number of observations, number of mismatches, empirical quality score).

#### input
- `bamfile:file`: A BAM file containing data that needs to be recalibrated.

#### output
- `outfile:file`: A GATKReport file with many tables:
- - The list of arguments
- - The quantized qualities table
- - The recalibration table by read group
- - The recalibration table by quality score
- - The recalibration table for all the optional covariates

#### args
- `bin`:     The gatk executable, default: "gatk -T BaseRecalibrator"
- `params`:  Other parameters for BaseRecalibrator, default: ""
- `reffile`: The reference file

#### requires
- [GATK](https://software.broadinstitute.org/gatk)


###  pPrintReads   
#### description
- PrintReads is a generic utility tool for manipulating sequencing data in SAM/BAM format. It can dynamically merge the contents of multiple input BAM files, resulting in merged output sorted in coordinate order. It can also optionally filter reads based on various read properties such as read group tags using the `--read_filter/-rf` command line argument (see documentation on read filters for more information).
- Note that when PrintReads is used as part of the Base Quality Score Recalibration workflow, it takes the `--BQSR` engine argument, which is listed under Inherited Arguments > CommandLineGATK below.

#### input
- `bamfile:file`: A BAM file.
- `recaltable:file`: The GATKReport file

#### output
- `outfile:file`: A single processed bam file.

#### args
- `bin`:     The gatk executable, default: "gatk -T PrintReads"
- `params`:  Other parameters for PrintReads, default: ""
- `reffile`: The reference file

#### requires
- [GATK](https://software.broadinstitute.org/gatk)


###  pHaplotypeCaller 
#### description
- PrintReads is a generic utility tool for manipulating sequencing data in SAM/BAM format. It can dynamically merge the contents of multiple input BAM files, resulting in merged output sorted in coordinate order. It can also optionally filter reads based on various read properties such as read group tags using the `--read_filter/-rf` command line argument (see documentation on read filters for more information).
- Note that when PrintReads is used as part of the Base Quality Score Recalibration workflow, it takes the `--BQSR` engine argument, which is listed under Inherited Arguments > CommandLineGATK below.

#### input
- `bamfile:file`: A BAM file.

#### output
- `outfile:file`: Either a VCF or gVCF file with raw, unfiltered SNP and indel calls.

#### args
- `bin`:     The gatk executable, default: "gatk -T HaplotypeCaller"
- `params`:  Other parameters for HaplotypeCaller, default: ""
- `reffile`: The reference file

#### requires
- [GATK](https://software.broadinstitute.org/gatk)


###  pSelectVariants
#### description
- Often, a VCF containing many samples and/or variants will need to be subset in order to facilitate certain analyses (e.g. comparing and contrasting cases vs. controls; extracting variant or non-variant loci that meet certain requirements, displaying just a few samples in a browser like IGV, etc.). SelectVariants can be used for this purpose.
- There are many different options for selecting subsets of variants from a larger callset:
- - Extract one or more samples from a callset based on either a complete sample name or a pattern match.
- - Specify criteria for inclusion that place thresholds on annotation values, e.g. "DP > 1000" (depth of coverage greater than 1000x), "AF < 0.25" (sites with allele frequency less than 0.25). These - criteria are written as "JEXL expressions", which are documented in the article about using JEXL expressions.
- - Provide concordance or discordance tracks in order to include or exclude variants that are also present in other given callsets.
- - Select variants based on criteria like their type (e.g. INDELs only), evidence of mendelian violation, filtering status, allelicity, and so on.
- There are also several options for recording the original values of certain annotations that are recalculated when a subsetting the new callset, trimming alleles, and so on.

#### input
- `vcffile:file`: A variant call set from which to select a subset.

#### output
- `outfile:file`: A new VCF file containing the selected subset of variants.

#### args
- `bin`:     The gatk executable, default: "gatk -T SelectVariants"
- `params`:  Other parameters for SelectVariants, default: ""
- `reffile`: The reference file

#### requires
- [GATK](https://software.broadinstitute.org/gatk)


###  pVariantFiltration
#### description
- This tool is designed for hard-filtering variant calls based on certain criteria. Records are hard-filtered by changing the value in the FILTER field to something other than PASS. Filtered records will be preserved in the output unless their removal is requested in the command line.
- The most common way of specifying filtering criteria is by using JEXL queries. See the article on JEXL expressions in the documentation Guide for detailed information and examples.

#### input
- `vcffile:file`: A variant set to filter.

#### output
- `outfile:file`: A filtered VCF.

#### args
- `bin`:     The gatk executable, default: "gatk -T VariantFiltration"
- `params`:  Other parameters for VariantFiltration, default: ""
- `reffile`: The reference file

#### requires
- [GATK](https://software.broadinstitute.org/gatk)


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
```
		             -(0.5 + 4*di/d0)
	PC = sum (pi * e                  )
```
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
```
		|--------- window ----------|
		|---- d0 -----|
		|--- 50K --- TSS --- 50K ---|
		     ^ (peak center)
		     |-- di --|
```


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
```
	pIntersectGMT.input = {pIntersectGMT.input: channel.create([(gmtfile1, gmtfile2)])}
	pIntersectGMT2 = pIntersectGMT.copy()
	pIntersectGMT2.depends = pIntersectGMT
	pIntersectGMT2.input   = {pIntersectGMT2.input.keys()[0]: lambda ch: ch.insert(0, gmtfile3)}
```

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
```
	pUnionGMT.input = {pIntersectGMT.input: channel.create([(gmtfile1, gmtfile2)])}
	pUnionGMT2 = pIntersectGMT.copy()
	pUnionGMT2.depends = pIntersectGMT
	pUnionGMT2.input   = {pUnionGMT.input.keys()[0]: lambda ch: ch.insert(0, gmtfile3)}
```

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
- - format: `<Probe name>\t<genotype>`
- - `<genotype>` = 0: AA, 1: AB, 2: BB

#### requires
- [bioconductor-crlmm](http://bioconductor.org/packages/release/bioc/html/crlmm.html)


###  pGenoToAvInput
#### description
- Convert the genotype called by pSNP6Genotype to [ANNOVAR input file](http://annovar.openbioinformatics.org/en/latest/user-guide/input/#annovar-input-file) using dbSNP identifiers.	

#### input
- `genofile:file`: the genofile generated by pSNP6Genotype, must be sorted by probe names
- `annofile:flie`: the annotation file downloaded from http://www.affymetrix.com/support/technical/annotationfilesmain.affx
- - Could be in .gz format

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
- `degfile:file`: the output file containing DEGs

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

