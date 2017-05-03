<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [Documentation for bioprocs v0.0.1](#documentation-for-bioprocs-v001)
  - [WXS](#wxs)
    - [pTrimmomaticPE](#ptrimmomaticpe)
    - [pTrimmomaticSE](#ptrimmomaticse)
    - [pAlignPEByBWA](#palignpebybwa)
    - [pAlignSEByBWA](#palignsebybwa)
    - [pAlignPEByNGM](#palignpebyngm)
    - [pAlignSEByNGM](#palignsebyngm)
    - [pCNVnator](#pcnvnator)
    - [pVcf2Stat](#pvcf2stat)
    - [pVcfStat2Mat](#pvcfstat2mat)
    - [pVcfStats](#pvcfstats)
    - [pCoverageByBamstats](#pcoveragebybamstats)
    - [pPlotBamstats](#pplotbamstats)
    - [pMutSig](#pmutsig)
    - [pSnpEff2Maf](#psnpeff2maf)
  - [CLUSTER](#cluster)
    - [pDist2Coords](#pdist2coords)
    - [pDecideK](#pdecidek)
    - [pKMeans](#pkmeans)
    - [pPamk](#ppamk)
    - [pClara](#pclara)
    - [pMClust](#pmclust)
    - [pAPCluster](#papcluster)
    - [pHClust](#phclust)
  - [PICARD](#picard)
    - [pMarkDuplicates](#pmarkduplicates)
    - [pAddOrReplaceReadGroups](#paddorreplacereadgroups)
    - [pCreateSequenceDictionary](#pcreatesequencedictionary)
    - [pCollectWgsMetrics](#pcollectwgsmetrics)
    - [pSortSam](#psortsam)
    - [pIndexBam](#pindexbam)
  - [WEB](#web)
    - [pDownloadPost](#pdownloadpost)
    - [pDownloadGet](#pdownloadget)
  - [TCGA](#tcga)
    - [pSample2SubmitterID](#psample2submitterid)
    - [pConvertExpFiles2Matrix](#pconvertexpfiles2matrix)
    - [pConvertMutFiles2Matrix](#pconvertmutfiles2matrix)
  - [ALGORITHM](#algorithm)
    - [pRWR](#prwr)
  - [GATK](#gatk)
    - [pRealignerTargetCreator](#prealignertargetcreator)
    - [pIndelRealigner](#pindelrealigner)
    - [pBaseRecalibrator](#pbaserecalibrator)
    - [pPrintReads](#pprintreads)
    - [pHaplotypeCaller](#phaplotypecaller)
    - [pSelectVariants](#pselectvariants)
    - [pVariantFiltration](#pvariantfiltration)
    - [pMuTect2](#pmutect2)
  - [GENOMEPLOT](#genomeplot)
    - [pGeneTrack](#pgenetrack)
    - [pAnnoTrack](#pannotrack)
    - [pDataTrack](#pdatatrack)
    - [pUcscTrack](#pucsctrack)
    - [pGenomePlot](#pgenomeplot)
  - [TFBS](#tfbs)
    - [pMotifScanByMEME](#pmotifscanbymeme)
    - [pMEMEMDB2Gene](#pmememdb2gene)
  - [SNPEFF](#snpeff)
    - [pAnn](#pann)
    - [pCount](#pcount)
  - [COMMON](#common)
    - [pSort](#psort)
  - [CHIPSEQ](#chipseq)
    - [pPeakToRegPotential](#ppeaktoregpotential)
  - [GSEA](#gsea)
    - [pMTarget2GTargetMat](#pmtarget2gtargetmat)
    - [pIntersectGMT](#pintersectgmt)
    - [pUnionGMT](#puniongmt)
    - [pSSGSEA](#pssgsea)
  - [SEQ](#seq)
    - [pConsv](#pconsv)
  - [SNPARRAY](#snparray)
    - [pSNP6Genotype](#psnp6genotype)
    - [pGenoToAvInput](#pgenotoavinput)
  - [DEG](#deg)
    - [pCallByLimmaFromMatrix](#pcallbylimmafrommatrix)
    - [pCallByLimmaFromFiles](#pcallbylimmafromfiles)
  - [NCBI](#ncbi)
  - [BEDTOOLS](#bedtools)
    - [pGetfasta](#pgetfasta)
    - [pClosest](#pclosest)
    - [pFlank](#pflank)
    - [pIntersect](#pintersect)
    - [pMakewindows](#pmakewindows)
    - [pMerge](#pmerge)
    - [pMultiinter](#pmultiinter)
    - [pRandom](#prandom)
    - [pShift](#pshift)
    - [pShuffle](#pshuffle)
    - [pSubtract](#psubtract)
    - [pWindow](#pwindow)
    - [pGenomecov](#pgenomecov)

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


###  pAlignPEByNGM
#### description
- Align paired-end reads to reference genome using NextGenMap

#### input
- `infile1:file`: read file 1 (fastq, or fastq gzipped)
- `infile2:file`: read file 2 (fastq, or fastq gzipped)

#### output
- `outfile:file`: The output sam/bam file

#### args
- `bin`:    The NextGenMap executable, default: ngm
- `params`: Other params for ngm, default: "--rg-id ngm --rg-sm sample"
- `nthread`: 1
- `reffile`: The reference file
- `outtype`: sam or bam, default: bam

#### requires
- [NextGenMap](https://github.com/Cibiv/NextGenMap/wiki)


###  pAlignSEByNGM
#### description
- Align single-end reads to reference genome using NextGenMap

#### input
- `infile:file`: read file (fastq, or fastq gzipped)

#### output
- `outfile:file`: The output sam/bam file

#### args
- `bin`:    The NextGenMap executable, default: ngm
- `params`: Other params for ngm, default: "--rg-id ngm --rg-sm sample"
- `nthread`: 1
- `reffile`: The reference file
- `outtype`: sam or bam, default: bam

#### requires
- [NextGenMap](https://github.com/Cibiv/NextGenMap/wiki)


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


###  pVcf2Stat
#### description
- Convert vcf to stat files for pVcfStats

#### input
- `vcffile:file`: The vcf file

#### output
- `outfile:file`: The stat file

#### args
- `chroms`: SNPs on chromosomes to consider, default: "*" (all chroms)
- - use "chr1-22, chrX, chrY" for chr1 to chr22, chrX and chrY

#### requires
- [`pyvcf`](https://github.com/jamescasbon/PyVCF)


###  pVcfStat2Mat
#### description
- Convert vcf stat files to matrix for clustering, rownames are sample names
- First 2 lines TiTv and Heterozygosity will be ignored.

#### input
- `statdir:file`: The input directory containing vcf stat files

#### output
- `outfile:file`: The matrix file


###  pVcfStats
#### description
- Calculate sample/snp call rate and heterozygosity, TiTv ratio from single-sample vcfs using result from pVcf2Stat

#### input
- `statdir:file`:   The directory containing vcf stat files

#### output
- `outsample:file`: The report of call rate for each sample
- `figsample:file`: The bar chat of sample call rates
- `outsnp:file`:    The report of call rate for each snp
- `figsnp:file`:    The bar chat of snp call rates
- `outhetero:file`: The report of heterozygosity of each sample
- `fighetero:file`: The bar chat of snp call rates
- `outtitv:file`:   The report of transition/transversion ratio of each sample
- `figtitv:file`:   The bar chat of transition/transversion ratio


###  pCoverageByBamstats
#### description
- Use `bamstats` to calculate coverage for bam file

#### input
- `infile:file`:  The bam file

#### output
- `outfile:file`:    The report of coverage for the bam file

#### args
- `bin`: The `bamstats` executable, default: "bamstats"
- `params`: Other parameters for `bamstats`, default: ""

#### requires
- [bamstats](http://bamstats.sourceforge.net/)


###  pPlotBamstats
#### description
- Plot coverage use output files generated by `bamstats` or `wxs.pCoverageByBamstats`

#### input
- `indir:file`: The directory containing bamstats output files

#### args
- `chroms`: Chromosomes to plot. Default: "" (all chroms)
- - Note: Whether to have "chr" prefix or not depends on your reference when mapping.
- - You can do a scope assignment: "chr1-chr22, chrX, chrY"

#### output
- `outdir:file`: The directory containing output figures


###  pMutSig
#### description
- MutSig stands for "Mutation Significance".  MutSig analyzes lists of mutations discovered in DNA sequencing, to identify genes that were mutated more often than expected by chance given background mutation processes.
- 
- For more information, see Lawrence, M. et al. Mutational heterogeneity in cancer and the search for new cancer-associated genes. Nature 499, 214-218 (2013).
- 
- See [dcumentation](http://archive.broadinstitute.org/cancer/cga/mutsig_run)

#### input
- `maffile:file`: mutation table
- `cvgfile:file`: coverage table
- `cvrfile:file`: covariates table
- `mutdict:file`: mutation_type_dictionary_file
- `chrdir:file`:  chr_files_hg18 or chr_files_hg19 

#### output
- `outfile:file`: The output file

#### args
- `bin`: The path to `run_MutSigCV.sh`, default: 'mutsig'
- `mcr`: The Matlab MCR path

#### requires
- [MutSing](http://archive.broadinstitute.org/cancer/cga/mutsig_download)


###  pSnpEff2Maf
#### description
- Convert a snpEff-annotated somatic mutation vcf file (with normal and tumor samples) to [maf](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification) file

#### input
- `vcffile:file`: snpEff-annotated vcf file

#### output
- `outfile:file`: The maf file

#### args
- `params`: A dict to specify a constant for columns. For example: `proc.args.params['NCBI_Build'] = 'hg19'` will set column 'NCBI_Build' as 'hg19' for all records. Default: {}
- - Keys could be one of these: ['Hugo_Symbol', 'Entrez_Gene_Id', 'Center', 'NCBI_Build', 'Chromosome', 'Start_Position', 'End_Position', 'Strand', 'Variant_Classification', 'Variant_Type', 'Reference_Allele', 'Tumor_Seq_Allele1', 'Tumor_Seq_Allele2', 'dbSNP_RS', 'dbSNP_Val_Status', 'Tumor_Sample_Barcode', 'Matched_Norm_Sample_Barcode', 'Match_Norm_Seq_Allele1', 'Match_Norm_Seq_Allele2', 'Tumor_Validation_Allele1', 'Tumor_Validation_Allele2', 'Match_Norm_Validation_Allele1', 'Match_Norm_Validation_Allele2', 'Verification_Status', 'Validation_Status', 'Mutation_Status', 'Sequencing_Phase', 'Sequence_Source', 'Validation_Method', 'Score', 'BAM_File', 'Tumor_Sample_UUID', 'Matched_Norm_Sample_UUID']
- `normal`: The normal sample position in `record.samples` from `pyvcf`, 0-based. Default: 0
- `tumor`:  The tumor sample position in `record.samples` from `pyvcf`, 0-based. Default: 1

#### requires
- [pyvcf](https://github.com/jamescasbon/PyVCF)


## CLUSTER

###  pDist2Coords
#### description
- Convert a distance matrix to 2D coordinates, using multidimensional scaling

#### input
- `infile:file`: The distance matrix, could be a full distance matrix, a triangle matrix or a pair-wise distance file
- - full dist matrix (full):
```
		s1	s2	s3
	s1	0	1	1
	s2	1	0	1
	s3	1	1	0
```
- - triangle matrix (triangle), could be also lower triangle
```
		s1	s2	s3
	s1	0	1	1
	s2		0	1
	s3			0
```
- - pair-wise (pair): (assuming auto-pair-wise distance = 0, that is: `s1	s1	0`)
```
	s1	s2	1
	s1	s3	1
	s2	s3	1
```
- - Both rownames and header of `full` and `triangle` can be omitted, just set `proc.args.rownames = "NULL"` and `proc.args.header = False`

#### output
- `outfile:file`: The output coordinate file

#### args
- `informat`: The format of the input file: full, triangle or pair. Default: full
- `rownames`: The `row.names` argument for `read.table`, default: 1
- `header`:   The `header` argument for `read.table` to read the input file, default: True.
- `k`:        How many dimension? Default: 2 (R^2)


###  pDecideK
#### description
- Decide number of clusters using different methods

#### input
- `infile:file`: the coordinates file, if all you have is a distance/similarity file, convert it to coordinates file using `pDist2Coords`

#### output
- `kfile:file`: the output file with `K`

#### args
- `method`:                         The method used to determine K
- - `elbow:<ev.thres>:<inc.thres>`: Look for a bend or elbow in the sum of squared error (SSE) scree plot, see [ref](https://artax.karlin.mff.cuni.cz/r-help/library/GMD/html/elbow.html). Default: `elbow` = `elbow:.95:01`
- - `pamk:<min>:<max>`:             You can do partitioning around medoids to estimate the number of clusters using the pamk function in the fpc package. Default: `pamk` = `pamk:2:15`
- - `calinski:<min>:<max>`:         Calinski criterion. Default: `calinski` means `calinski:2:15`
- - `mclust:<min>:<max>`:           Determine the optimal model and number of clusters according to the Bayesian Information Criterion for expectation-maximization, initialized by hierarchical clustering for parameterized Gaussian mixture models. [Ref1](http://www.stat.washington.edu/research/reports/2006/tr504.pdf
- #), [Ref2](http://www.jstatsoft.org/v18/i06/paper). Default: `mclust` = `mclust:2:15`
- - `ap`:                           Affinity propagation (AP) clustering, see [ref](http://dx.doi.org/10.1126/science.1136800)
- - `gap:<min>:<max>`:              Gap Statistic for Estimating the Number of Clusters. Default: `gap:2:10`
- - `nbclust`:                      The [NbClust package](http://cran.r-project.org/web/packages/NbClust/index.html) provides 30 indices to determine the number of clusters in a dataset.
- `rownames`:                       The `row.names` for `read.table` to read the input file, default: 1.
- `header`:                         The `header` argument for `read.table` to read the input file, default: True.
- `seed`:                           The seed for randomization, default: 0.

#### requires
- [`r-cluster`](https://cran.r-project.org/web/packages/cluster/index.html) if `gap` method used
- [`r-GMD`](https://cran.r-project.org/web/packages/GMD/index.html) if `elbow` method userd
- [`r-fpc`](https://cran.r-project.org/web/packages/fpc/index.html) if `pamk` method used
- [`r-vegan`](https://cran.r-project.org/web/packages/vegan/index.html) if `calinski` method used
- [`r-mclust`](https://cran.r-project.org/web/packages/mclust/index.html) if `mclust` method used
- [`r-apcluster`](https://cran.r-project.org/web/packages/apcluster/index.html) if `ap` method used
- [`r-NbClust`](https://cran.r-project.org/web/packages/NbClust/index.html) if `nbclust` method used


###  pKMeans
#### description
- Do k-means clustering

#### input
- `infile:file`:    The input coordinates of the points.
- `k`:              Number of clusters, it could also be a file with the number in it.

#### output
- `outdir:dir`: The output of final results

#### args
- `rownames`:       The `row.names` for `read.table` to read the input file, default: 1.
- `header`:         The `header` argument for `read.table` to read the input file, default: True.
- `algorithm`:      The `algorithm` argument for `kmeans`, default "Hartigan-Wong" (could also be "Lloyd", "Forgy", "MacQueen")
- `niter`:          The `max.iter` argument for `kmeans`, default: 10.
- `nstart`:         The `nstart` argument for `kmeans`, default: 25.
- `caption`:        The caption for the `fviz_cluster`, default: "K-means with K=%k%".

#### requires
- [`r-factoextra`](https://cran.r-project.org/web/packages/factoextra/index.html)


###  pPamk
#### description
- Do clustering using [fpc::pamk](https://www.rdocumentation.org/packages/fpc/versions/2.1-10/topics/pamk)

#### input
- `infile:file`:  The input coordinate file

#### output
- `outdir:dir`:   The output directory

#### args
- `rownames`:     The `row.names` for `read.table` to read the input file, default: 1.
- `header`:       The `header` argument for `read.table` to read the input file, default: True.
- `min`:          The min # clusters to try, default: 2
- `max`:          The max # clusters to try, default: 15
- `caption`:      The caption for the `fviz_cluster`, default: "Partitioning Around Medoids (K=%K%)".
- `seed`:         The seed for randomization, default: 0.

#### requires
- [`r-factoextra`](https://cran.r-project.org/web/packages/factoextra/index.html)
- [`r-fpc`](https://cran.r-project.org/web/packages/fpc/index.html)


###  pClara
#### description
- CLARA is a partitioning method used to deal with much larger data sets (more than several thousand observations) in order to reduce computing time and RAM storage problem.

#### input
- `infile:file`:  The input coordinate file
- `k`:            Number of clusters, it could also be a file with the number in it.

#### output
- `outdir:dir`:   The output of final results

#### args
- `rownames`:     The `row.names` for `read.table` to read the input file, default: 1.
- `header`:       The `header` argument for `read.table` to read the input file, default: True.
- `samples`:      The `samples` argument for `clara`, default: 5.
- `caption`:      The caption for the `fviz_cluster`, default: "CLARA Clustering with K=%k%".

#### requires
- [`r-cluster`](https://cran.r-project.org/web/packages/cluster/index.html)
- [`r-factoextra`](https://cran.r-project.org/web/packages/factoextra/index.html)


###  pMClust
#### description
- Use `r-mclust` to do clustering. Current just do simple clustering with the package

#### input
- `infile:file`:  The input a coordinate file

#### output
- `outdir:dir`:   The output of final results

#### args
- `rownames`:     The `row.names` for `read.table` to read the input file, default: 1.
- `header`:       The `header` argument for `read.table` to read the input file, default: True.
- `caption`:      The caption for the `fviz_cluster`, default: "CLARA Clustering".
- `min`:          The min # clusters to try, default: 2
- `max`:          The max # clusters to try, default: 15

#### requires
- [`r-mclust`](https://cran.r-project.org/web/packages/mclust/index.html)
- [`r-factoextra`](https://cran.r-project.org/web/packages/factoextra/index.html)


###  pAPCluster
#### description
- Use `r-apcluster` to do clustering. 

#### input
- `infile:file`:  The input a coordinate file

#### output
- `outdir:dir`:   The output of final results

#### args
- `rownames`:     The `row.names` for `read.table` to read the input file, default: 1.
- `header`:       The `header` argument for `read.table` to read the input file, default: True.
- `caption`:      The caption for the `fviz_cluster`, default: "APClustering".

#### requires
- [`r-apcluster`](https://cran.r-project.org/web/packages/apcluster/index.html)
- [`r-factoextra`](https://cran.r-project.org/web/packages/factoextra/index.html)


###  pHClust
#### description
- Do hierarchical clustering.

#### input
- `infile:file`: The input files with variants as rows, features as columns.
- - NOTE: clustering is performed on rows, rownames are the leaf labels.

#### output
- `outdir:dir`:  The result directory, containing:
- - `hclust.merge.txt`: including merge and height information
- - `hclust.order.txt`: including order and labels information
- - `hclust.png`:       the dendrogram plot

#### args
- `fast`:     whether to use `fastcluster` package or not, default: False
- `gg`:       whether to use `ggdendro` or not, default: False
- `rownames`: The `row.names` for `read.table` to read the input file, default: 1.
- `header`:   The `header` argument for `read.table` to read the input file, default: True.
- `method`:   Which method to use for `hclust`. Default: "complete" (use `?hclust` to check all availables)
- `rotate`:   Which to rotate the plot or not. Default: False

#### requires
- [`r-fastcluster`](https://cran.r-project.org/web/packages/fastcluster/index.html) if `proc.args.fast` is True
- [`r-ggdendro`](https://cran.r-project.org/web/packages/ggdendro/index.html) if `proc.args.gg` is True


## PICARD

###  pMarkDuplicates
#### description
- Identifies duplicate reads.
- This tool locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating from a single fragment of DNA. Duplicates can arise during sample preparation e.g. library construction using PCR. See also EstimateLibraryComplexity for additional notes on PCR duplication artifacts. Duplicate reads can also result from a single amplification cluster, incorrectly detected as multiple clusters by the optical sensor of the sequencing instrument. These duplication artifacts are referred to as optical duplicates.
- The MarkDuplicates tool works by comparing sequences in the 5 prime positions of both reads and read-pairs in a SAM/BAM file. An BARCODE_TAG option is available to facilitate duplicate marking using molecular barcodes. After duplicate reads are collected, the tool differentiates the primary and duplicate reads using an algorithm that ranks reads by the sums of their base-quality scores (default method).
- The tool's main output is a new SAM or BAM file, in which duplicates have been identified in the SAM flags field for each read. Duplicates are marked with the hexadecimal value of 0x0400, which corresponds to a decimal value of 1024. If you are not familiar with this type of annotation, please see the following [blog post](https://www.broadinstitute.org/gatk/blog?id=7019) for additional information.
- Although the bitwise flag annotation indicates whether a read was marked as a duplicate, it does not identify the type of duplicate. To do this, a new tag called the duplicate type (DT) tag was recently added as an optional output in the 'optional field' section of a SAM/BAM file. Invoking the TAGGING_POLICY option, you can instruct the program to mark all the duplicates (All), only the optical duplicates (OpticalOnly), or no duplicates (DontTag). The records within the output of a SAM/BAM file will have values for the 'DT' tag (depending on the invoked TAGGING_POLICY), as either library/PCR-generated duplicates (LB), or sequencing-platform artifact duplicates (SQ). This tool uses the READ_NAME_REGEX and the OPTICAL_DUPLICATE_PIXEL_DISTANCE options as the primary methods to identify and differentiate duplicate types. Set READ_NAME_REGEX to null to skip optical duplicate detection, e.g. for RNA-seq or other data where duplicate sets are extremely large and estimating library complexity is not an aim. Note that without optical duplicate counts, library size estimation will be inaccurate.
- MarkDuplicates also produces a metrics file indicating the numbers of duplicates for both single- and paired-end reads.
- The program can take either coordinate-sorted or query-sorted inputs, however the behavior is slightly different. When the input is coordinate-sorted, unmapped mates of mapped records and supplementary/secondary alignments are not marked as duplicates. However, when the input is query-sorted (actually query-grouped), then unmapped mates and secondary/supplementary reads are not excluded from the duplication test and can be marked as duplicate reads.
- If desired, duplicates can be removed using the REMOVE_DUPLICATE and REMOVE_SEQUENCING_DUPLICATES options.

#### input
- `infile:file`:  The bam file 

#### output
- `outfile:file`: The marked bam file

#### args
- `bin`:     The picard executable, default: "picard MarkDuplicates"
- `params`:  Other parameters for picard MarkDuplicates, default: ""

#### requires
- [picard](https://broadinstitute.github.io/picard/)


###  pAddOrReplaceReadGroups
#### description
- Replace read groups in a BAM file.This tool enables the user to replace all read groups in the INPUT file with a single new read group and assign all reads to this read group in the OUTPUT BAM file.
- For more information about read groups, see the [GATK Dictionary entry](https://www.broadinstitute.org/gatk/guide/article?id=6472). 
- This tool accepts INPUT BAM and SAM files or URLs from the Global Alliance for Genomics and Health (GA4GH) (see http://ga4gh.org/#/documentation).

#### input
- `infile:file`:  The bam file 

#### output
- `outfile:file`: The bam file with read group added

#### args
- `bin`:     The picard executable, default: "picard AddOrReplaceReadGroups"
- `params`:  Other parameters for picard AddOrReplaceReadGroups, default: "RGID=4 RGLB=lib1 RGPL=illumina  RGPU=unit1 RGSM=20"

#### requires
- [picard](https://broadinstitute.github.io/picard/)


###  pCreateSequenceDictionary
#### description
- Creates a sequence dictionary for a reference sequence. This tool creates a sequence dictionary file (with ".dict" extension) from a reference sequence provided in FASTA format, which is required by many processing and analysis tools. The output file contains a header but no SAMRecords, and the header contains only sequence records.
- The reference sequence can be gzipped (both .fasta and .fasta.gz are supported).

#### input
- `infile:file`:  The fasta file 

#### output
- `outfile:file`: The same fasta file, but with dict file created

#### args
- `bin`:     The picard executable, default: "picard CreateSequenceDictionary"
- `params`:  Other parameters for picard CreateSequenceDictionary, default: ""

#### requires
- [picard](https://broadinstitute.github.io/picard/)


###  pCollectWgsMetrics
#### description
- Collect metrics about coverage and performance of whole genome sequencing (WGS) experiments.
- This tool collects metrics about the fractions of reads that pass base- and mapping-quality filters as well as coverage (read-depth) levels for WGS analyses. Both minimum base- and mapping-quality values as well as the maximum read depths (coverage cap) are user defined.
- Note: Metrics labeled as percentages are actually expressed as fractions!

#### input
- `infile:file`:  The bam file 

#### output
- `outfile:file`: The metrics file

#### args
- `bin`:     The picard executable, default: "picard CollectWgsMetrics"
- `params`:  Other parameters for `picard CollectWgsMetrics`, default: ""
- `reffile`: The reference file, default: ""

#### requires
- [picard](https://broadinstitute.github.io/picard/)


###  pSortSam
#### description
- Use `picard SortSam` to sort sam or bam file

#### input
- `infile:file`:  The sam or bam file to be sorted

#### output
- `outfile:file`: The sorted sam or bam file

#### args
- `bin`:     The picard executable, default: "picard SortSam"
- `order`:   The sort order, default: coordinate. Possible: unsorted, queryname, coordinate, duplicate
- `outtype`: The type of output file, sam or bam. Default: bam
- `params`:  Other parameters for `picard SortSame`, default: "-Xms1g -Xmx8g"

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
- `params`:  Other parameters for picard , default: "-Xms1g -Xmx8g"

#### requires
- [picard](http://broadinstitute.github.io/picard/command-line-overview.html)


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
- `bin-samtools`: The samtools executable, default: samtools

#### requires
- [GATK](https://software.broadinstitute.org/gatk)
- [samtools](http://www.htslib.org/) if `reffile` is not indexed or `infile` is not indexed.


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
- `bin-samtools`: The samtools executable, default: samtools

#### requires
- [GATK](https://software.broadinstitute.org/gatk)
- [samtools](http://www.htslib.org/) if `reffile` is not indexed or `infile` is not indexed.


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
- `reffile`: The reference file, required
- `knownSites`: The known polymorphic sites to mask out, required

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


###  pMuTect2
#### description
- MuTect2 is a somatic SNP and indel caller that combines the DREAM challenge-winning somatic genotyping engine of the original MuTect ([Cibulskis et al., 2013](http://www.nature.com/nbt/journal/v31/n3/full/nbt.2514.html)) with the assembly-based machinery of HaplotypeCaller. The basic operation of MuTect2 proceeds similarly to that of the HaplotypeCaller.
- NOTE: only Tumor/Normal variant calling implemented in bioprocs

#### input
- `tumor:file`:  the tumor bam file
- `normal:file`: the normal bam file

#### output
- `outfile:file`: The vcf file containing somatic mutations

#### args
- `bin`:     The gatk executable, default: "gatk -T MuTect2"
- `params`:  Other parameters for MuTect2, default: ""
- `reffile`: The reference file
- `bin-samtools`: the samtools executable, default: samtools

#### requires
- [GATK](https://software.broadinstitute.org/gatk)
- [samtools](http://www.htslib.org/) if `reffile` is not indexed or `infile` is not indexed.


## GENOMEPLOT

###      pGeneTrack
#### description
- Generate the gene track using ucsc data source

#### input
- `name`:   The name of the track
- `genome`: The genome
- `chrom`:  The chromosome
- `from`:   The start
- `to`:     The end

#### output
- `outfile:file`: The file to save the track

#### args
- use `displayPars(UcscTrack(genome="mm9", chromosome="chrM", track="knownGene"))` to see all available args

#### requires
- [r-Gviz](https://rdrr.io/bioc/Gviz)


###      pAnnoTrack
#### description
- Generate annotation track

#### input
- `infile:file`: the file for the track
- `name`:        the name of the track
- `genome`:      the genome
- `chrom`:       the chromosome
- `from`:        the start position to display
- `to`:          the end position to display

#### output
- `outfile:file`:the dumped track

#### args
- use `displayPars(AnnotationTrack())` to see all available args.

#### requires
- [r-Gviz](https://rdrr.io/bioc/Gviz)


###      pDataTrack
#### description
- The data track of Gviz

#### input
- `infile:file`: the file for the track
- `name`:        the name of the track
- `genome`:      the genome
- `chrom`:       the chromosome
- `from`:        the start position to display
- `to`:          the end position to display

#### output
- `outfile:file`:the rds file for the track
- `gout`:        the genome
- `cout`:        the chromosome
- `fout`:        the start
- `tout`:        the end

#### args
- See `displayPars(DataTrack())` for all available display params
- Quote all params!

#### requires
- [r-Gviz](https://rdrr.io/bioc/Gviz/man/DataTrack-class.html)


###      pUcscTrack
#### description
- Generate track from ucsc

#### input
- `ucscTrack`:   the track to fetch from ucsc. [Avialable tracks](http://genome.ucsc.edu/cgi-bin/hgTables?command=start)
- `table`:       the table from ucsc. [Available table](http://genome.ucsc.edu/cgi-bin/hgTables?command=start)
- `gvizTrack`:   the object track to generate. One of "AnnotationTrack", "GeneRegionTrack", "DataTrack", "GenomeAxisTrack"
- `name`:        the name of the track
- `genome`:      the genome
- `chrom`:       the chromosome
- `from`:        the start position to display
- `to`:          the end position to display

#### output
- `outfile:file`:the dumped track

#### args
- use `displayPars(UcscTrack(genome="mm9", chromosome="chrM", track="knownGene"))` to see all available args.

#### requires
- [r-Gviz](https://rdrr.io/bioc/Gviz)


###      pGenomePlot
#### description
- plot the genomic features

#### input
- `trkfiles:files`: the list of track dumped files
- `genome`:         the genome
- `chrom`:          the chromosome
- `from`:           the start position to display
- `to`:             the end position to display

#### output
- `outfile:file`:   the figure

#### requires
- [r-Gviz](https://rdrr.io/bioc/Gviz)


## TFBS

###  pMotifScanByMEME
#### input
- `mfile:file`: The motif file
- `sfile:file`: The sequence file

#### output
- `outdir:file`: The output dir

#### args
- `bin-fimo`: The path of `fimo` executable, default: "fimo"
- `params`:   Other parameters for `fimo`

#### requires
- [`fimo` from MEME Suite](http://meme-suite.org/tools/fimo)


###  pMEMEMDB2Gene
#### input
- `memefile:file`: The meme motif downloaded from MEME website
- `species`:       The species, can be tax id, or common name
- - only can be one of fruitfly, mouse, human, frog, zebrafish, thale-cress, pig, rat, nematode
- - check them from http://http://mygene.info/v3/metadata

#### output
- `outfile:file`:  The output file containing the name pairs

#### requires
- [python-mygene](https://pypi.python.org/pypi/mygene)


## SNPEFF

###  pAnn
#### description
- This is the default command. It is used for annotating variant filed (e.g. VCF files).

#### input
- `infile:file`:  The input file 

#### output
- `outdir:file`: The directory containing output anntated file, snpEff_genes.txt and snpEff_summary.html

#### args
- `bin`:       The snpEff executable, default: "snpEff"
- `params`:    Other parameters for `snpEff`, default: "-Xms1g -Xmx4g -v"
- `genome`:    The genome used for annotation, default: "hg19"
- `informat`:  The format of input file [vcf or bed], default: "vcf"
- `outformat`: The format of output file [vcf, gatk, bed, bedAnn], default: "vcf"
- `csvStats`:  Whether to generate csv stats file, default: "{{infile.fn}}.stats.csv", set False to disable.
- `stats`:     The name of the summary file, default: "{{infile.fn}}.summary.html", set False to disable.

#### requires
- [snpEff](http://snpeff.sourceforge.net/SnpEff_manual.html)


###  pCount
#### description
- Count how many intervals (from a BAM, BED or VCF file) overlap with each genomic interval.

#### input
- `infiles:files`:  The input files

#### output
- `outdir:file`: The directory containg summary (html) file and output file

#### args
- `bin`:     The snpEff executable, default: "snpEff"
- `params`:  Other parameters for `snpEff`, default: "-Xms1g -Xmx4g -v"
- `genome`:  The genome used for counting, default: "hg19"

#### requires
- [snpEff](http://snpeff.sourceforge.net/SnpEff_manual.html)


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


## SEQ

###  pConsv
#### description
- Get the conservation scores of regions.
- It uses wigFix to find the conservation scores.
- But first you have to convert those wigFix.gz files to bigWig files using ucsc-wigToBigWig

#### input
- `bedfile:file`: The bedfile with regions in the same chromosome
- `bwdir:file`:   The bigwig directory, the bigwig files must be named as "chrN.*"
- - For example: `chr1.phyloP30way.bw`

#### output
- `outdir:dir`    The output file, containing the output file, random regions and their scores if pvalue is calulated

#### args
- `bin-bwtool`:   The path of bwtool executable. Default: `bwtool`
- `bin-bedtools`: The path of bedtools executable. Default: `bedtools`
- `calcp`:        Whether to calculate the pvalue for not for each region.
- `nperm`:        The number of permutations of pvalue is calcuated.
- `seed`:         The seed to generate the random regions.
- `gsize`:        The genome size file for generating random regions.
- - Can be a local file or a file hosted at UCSC: http://hgdownload.cse.ucsc.edu/goldenPath/<genome>/bigZips/<genome>.chrom.sizes

#### requires
- [bwtool](https://github.com/CRG-Barcelona/bwtool)
- [bedtools](http://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html) if `calcp` is `True`


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


## NCBI

## BEDTOOLS

###  pGetfasta
#### description
- `bedtools getfasta` extracts sequences from a FASTA file for each of the intervals defined in a BED/GFF/VCF file.

#### input
- `infile:file`: The input bed file
- `fafile:file`: The input fasta file

#### brings
- `fafile`: "{{fafile.fn}}.fa*i", The fasta index file

#### output
- `outfile:file`: The generated fasta file

#### args
- `bin`:     The bedtools executable, default: "bedtools"
- `params`:  Other parameters for `bedtools getfasta`, default: ""

#### requires
- [bedtools](http://bedtools.readthedocs.io/en/latest/index.html)


###  pClosest
#### description
- Similar to intersect, closest searches for overlapping features in A and B. In the event that no feature in B overlaps the current feature in A, closest will report the nearest (that is, least genomic distance from the start or end of A) feature in B. For example, one might want to find which is the closest gene to a significant GWAS polymorphism. Note that closest will report an overlapping feature as the closest that is, it does not restrict to closest non-overlapping feature. The following iconic cheatsheet summarizes the funcitonality available through the various optyions provided by the closest tool.

#### input
- `afile:file`:   The -a file
- `bfiles:files`: The -b files

#### output
- `outfile:file`: The result file

#### args
- `bin`:     The bedtools executable, default: "bedtools"
- `params`:  Other parameters for `bedtools closest`, default: ""

#### requires
- [bedtools](http://bedtools.readthedocs.io/en/latest/index.html)


###  pFlank
#### description
- `bedtools flank` will create two new flanking intervals for each interval in a BED/GFF/VCF file. Note that flank will restrict the created flanking intervals to the size of the chromosome (i.e. no start < 0 and no end > chromosome size).

#### input
- `infile:file`:  The input file
- `gfile:file`:   The genome size file

#### output
- `outfile:file`: The result file

#### args
- `bin`:     The bedtools executable, default: "bedtools"
- `params`:  Other parameters for `bedtools flank`, default: ""

#### requires
- [bedtools](http://bedtools.readthedocs.io/en/latest/index.html)


###  pIntersect
#### description
- By far, the most common question asked of two sets of genomic features is whether or not any of the features in the two sets overlap with one another. This is known as feature intersection. bedtools intersect allows one to screen for overlaps between two sets of genomic features. Moreover, it allows one to have fine control as to how the intersections are reported. bedtools intersect works with both BED/GFF/VCF and BAM files as input.

#### input
- `afile:file`:   The a file
- `bfiles:files`: The b files

#### output
- `outfile:file`: The result file

#### args
- `bin`:     The bedtools executable, default: "bedtools"
- `params`:  Other parameters for `bedtools intersect`, default: ""

#### requires
- [bedtools](http://bedtools.readthedocs.io/en/latest/index.html)


###  pMakewindows
#### description
- Makes adjacent or sliding windows across a genome or BED file.

#### input
- `infile:file`: The input file

#### output
- `outfile:file`: The result file

#### args
- `bin`:     The bedtools executable, default: "bedtools"
- `informat`:The format of input file, whether is a "bed" file or "genome" size file. Default: "bed"
- `params`:  Other parameters for `bedtools makewindows`, default: ""

#### requires
- [bedtools](http://bedtools.readthedocs.io/en/latest/index.html)


###  pMerge
#### description
- `bedtools merge` combines overlapping or book-ended features in an interval file into a single feature which spans all of the combined features.

#### input
- `infile:file`: The input file

#### output
- `outfile:file`: The result file

#### args
- `bin`:     The bedtools executable, default: "bedtools"
- `params`:  Other parameters for `bedtools merge`, default: ""

#### requires
- [bedtools](http://bedtools.readthedocs.io/en/latest/index.html)


###  pMultiinter
#### description
- Identifies common intervals among multiple BED/GFF/VCF files.

#### input
- `infiles:files`: The input files

#### output
- `outfile:file`: The result file

#### args
- `bin`:     The bedtools executable, default: "bedtools"
- `params`:  Other parameters for `bedtools multiinter`, default: ""

#### requires
- [bedtools](http://bedtools.readthedocs.io/en/latest/index.html)


###  pRandom
#### description
- `bedtools random` will generate a random set of intervals in BED6 format. One can specify both the number (-n) and the size (-l) of the intervals that should be generated.

#### input
- `gfile:file`: The genome size file

#### output
- `outfile:file`: The result file

#### args
- `bin`:     The bedtools executable, default: "bedtools"
- `params`:  Other parameters for `bedtools random`, default: ""

#### requires
- [bedtools](http://bedtools.readthedocs.io/en/latest/index.html)


###  pShift
#### description
- `bedtools shift` will move each feature in a feature file by a user-defined number of bases. While something like this could be done with an awk '{OFS="\t" print $1,$2+<shift>,$3+<shift>}', bedtools shift will restrict the resizing to the size of the chromosome (i.e. no features before 0 or past the chromosome end).

#### input
- `infile:file`: The input file
- `gfile:file`:  The genome size file

#### output
- `outfile:file`: The result file

#### args
- `bin`:     The bedtools executable, default: "bedtools"
- `params`:  Other parameters for `bedtools shift`, default: ""

#### requires
- [bedtools](http://bedtools.readthedocs.io/en/latest/index.html)


###  pShuffle
#### description
- `bedtools shuffle` will randomly permute the genomic locations of a feature file among a genome defined in a genome file. One can also provide an exclusions BED/GFF/VCF file that lists regions where you do not want the permuted features to be placed. For example, one might want to prevent features from being placed in known genome gaps. shuffle is useful as a null basis against which to test the significance of associations of one feature with another.

#### input
- `infile:file`: The input file
- `gfile:file`:  The genome size file

#### output
- `outfile:file`: The result file

#### args
- `bin`:     The bedtools executable, default: "bedtools"
- `params`:  Other parameters for `bedtools shuffle`, default: ""

#### requires
- [bedtools](http://bedtools.readthedocs.io/en/latest/index.html)


###  pSubtract
#### description
- `bedtools subtract` searches for features in B that overlap A. If an overlapping feature is found in B, the overlapping portion is removed from A and the remaining portion of A is reported. If a feature in B overlaps all of a feature in A, the A feature will not be reported.

#### input
- `afile:file`: The a file
- `bfile:file`: The b file

#### output
- `outfile:file`: The result file

#### args
- `bin`:     The bedtools executable, default: "bedtools"
- `params`:  Other parameters for `bedtools subtract`, default: ""

#### requires
- [bedtools](http://bedtools.readthedocs.io/en/latest/index.html)


###  pWindow
#### description
- Similar to `bedtools intersect`, `window` searches for overlapping features in A and B. However, window adds a specified number (1000, by default) of base pairs upstream and downstream of each feature in A. In effect, this allows features in B that are near features in A to be detected.

#### input
- `afile:file`: The a file
- `bfile:file`: The b file

#### output
- `outfile:file`: The result file

#### args
- `bin`:     The bedtools executable, default: "bedtools"
- `params`:  Other parameters for `bedtools window`, default: ""

#### requires
- [bedtools](http://bedtools.readthedocs.io/en/latest/index.html)


###  pGenomecov
#### description
- `bedtools genomecov` computes histograms (default), per-base reports (-d) and BEDGRAPH (-bg) summaries of feature coverage (e.g., aligned sequences) for a given genome.
- NOTE: only bam file input implemented here.

#### input
- `infile:file`: The bam file

#### output
- `outfile:file`: The result file

#### args
- `bin`:     The bedtools executable, default: "bedtools"
- `params`:  Other parameters for `bedtools genomecov`, default: "-bg"

#### requires
- [bedtools](http://bedtools.readthedocs.io/en/latest/index.html)

