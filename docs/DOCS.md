<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [Documentation for bioprocs v0.0.1a](#documentation-for-bioprocs-v001a)
  - [CLUSTER](#cluster)
    - [pDist2Coords](#pdist2coords)
    - [pCluster](#pcluster)
    - [pMClust](#pmclust)
    - [pAPCluster](#papcluster)
    - [pHClust](#phclust)
  - [FASTX](#fastx)
    - [pFastqPESim](#pfastqpesim)
    - [pFastQC](#pfastqc)
    - [pFastMC](#pfastmc)
    - [pFastqPETrim](#pfastqpetrim)
    - [pFastqSETrim](#pfastqsetrim)
    - [pFastqSE2Sam](#pfastqse2sam)
    - [pFastqPE2Sam](#pfastqpe2sam)
  - [WXSDOWN](#wxsdown)
    - [pMutSig](#pmutsig)
    - [pVcf2Maf](#pvcf2maf)
    - [pMergeMafs](#pmergemafs)
    - [pMutsig4Plot](#pmutsig4plot)
    - [pMutPlot](#pmutplot)
    - [pCepip](#pcepip)
  - [VCF](#vcf)
    - [pVcfFilter](#pvcffilter)
    - [pVcfAnno](#pvcfanno)
    - [pCallRate](#pcallrate)
  - [WXSCALL](#wxscall)
    - [pCNVnator](#pcnvnator)
  - [GSEA](#gsea)
    - [pMTarget2GTargetMat](#pmtarget2gtargetmat)
    - [pIntersectGMT](#pintersectgmt)
    - [pUnionGMT](#puniongmt)
    - [pExpmat2Gct](#pexpmat2gct)
    - [pSampleinfo2Cls](#psampleinfo2cls)
    - [pSSGSEA](#pssgsea)
    - [pGSEA](#pgsea)
    - [pEnrichr](#penrichr)
    - [pTargetEnrichr](#ptargetenrichr)
  - [ALGORITHM](#algorithm)
    - [pRWR](#prwr)
  - [SNPINFO](#snpinfo)
    - [pSNP2Bed](#psnp2bed)
    - [pSNP2Avinput](#psnp2avinput)
  - [MARRAY](#marray)
    - [pCeldir2Matrix](#pceldir2matrix)
    - [aCelPat2Deg](#acelpat2deg)
    - [aCelPat2DegGSEA](#acelpat2deggsea)
  - [CHIPSEQ](#chipseq)
    - [pPeakToRegPotential](#ppeaktoregpotential)
  - [PLOT](#plot)
    - [pBoxplot](#pboxplot)
    - [pScatterPlot](#pscatterplot)
    - [pVenn](#pvenn)
  - [WXSANNO](#wxsanno)
    - [pSnpEff](#psnpeff)
  - [RNASEQ](#rnaseq)
    - [pExpdir2Matrix](#pexpdir2matrix)
    - [pBatchEffect](#pbatcheffect)
    - [pRawCounts2](#prawcounts2)
    - [pDeg](#pdeg)
  - [WEB](#web)
    - [pDownloadPost](#pdownloadpost)
    - [pDownloadGet](#pdownloadget)
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
  - [RANK](#rank)
    - [pRankProduct](#prankproduct)
  - [WXSPREP](#wxsprep)
    - [pTrimmomaticPE](#ptrimmomaticpe)
    - [pTrimmomaticSE](#ptrimmomaticse)
    - [pAlignPEByBWA](#palignpebybwa)
    - [pAlignSEByBWA](#palignsebybwa)
    - [pAlignPEByNGM](#palignpebyngm)
    - [pAlignSEByNGM](#palignsebyngm)
    - [pMergeBams](#pmergebams)
  - [DEG](#deg)
    - [pExpdirMatrix](#pexpdirmatrix)
    - [pDEGByEdgeR](#pdegbyedger)
    - [pMArrayLimma](#pmarraylimma)
    - [pRawCounts2](#prawcounts2-1)
  - [VCFNEXT](#vcfnext)
    - [pStats2Matrix](#pstats2matrix)
    - [pPlotStats](#pplotstats)
  - [PICARD](#picard)
    - [pMarkDuplicates](#pmarkduplicates)
    - [pAddOrReplaceReadGroups](#paddorreplacereadgroups)
    - [pCreateSequenceDictionary](#pcreatesequencedictionary)
    - [pCollectWgsMetrics](#pcollectwgsmetrics)
    - [pSortSam](#psortsam)
    - [pIndexBam](#pindexbam)
  - [GENOMEPLOT](#genomeplot)
    - [pGeneTrack](#pgenetrack)
    - [pAnnoTrack](#pannotrack)
    - [pDataTrack](#pdatatrack)
    - [pUcscTrack](#pucsctrack)
    - [pGenomePlot](#pgenomeplot)
  - [BED](#bed)
    - [pBedSort](#pbedsort)
    - [pBedIntersect](#pbedintersect)
    - [pBedCluster](#pbedcluster)
  - [RESOURCE](#resource)
    - [pTxt](#ptxt)
  - [PCA](#pca)
    - [pPCA](#ppca)
    - [pSelectPCs](#pselectpcs)
  - [TFBS](#tfbs)
    - [pMotifScanByMEME](#pmotifscanbymeme)
    - [pMEMEMDB2Gene](#pmememdb2gene)
  - [SAMBAM](#sambam)
    - [pSam2Bam](#psam2bam)
    - [pBamMarkdup](#pbammarkdup)
    - [pBamRecal](#pbamrecal)
    - [pBamReadGroup](#pbamreadgroup)
    - [pBamReorder](#pbamreorder)
    - [pBamMerge](#pbammerge)
    - [pBam2Gmut](#pbam2gmut)
    - [pBam2Cnv](#pbam2cnv)
    - [pBam2FastqPE](#pbam2fastqpe)
    - [pBam2FastqSE](#pbam2fastqse)
    - [pBam2Counts](#pbam2counts)
  - [GATK](#gatk)
    - [pRealignerTargetCreator](#prealignertargetcreator)
    - [pIndelRealigner](#pindelrealigner)
    - [pBaseRecalibrator](#pbaserecalibrator)
    - [pPrintReads](#pprintreads)
    - [pHaplotypeCaller](#phaplotypecaller)
    - [pSelectVariants](#pselectvariants)
    - [pVariantFiltration](#pvariantfiltration)
    - [pMuTect2](#pmutect2)
    - [pMuTect2Interval](#pmutect2interval)
  - [UTILS](#utils)
  - [WXS](#wxs)
  - [COMMON](#common)
    - [pSort](#psort)
    - [pFiles2Dir](#pfiles2dir)
    - [pFiles2List](#pfiles2list)
    - [pPat2Dir](#ppat2dir)
    - [pMergeFiles](#pmergefiles)
    - [pCbindList](#pcbindlist)
    - [pFile2Proc](#pfile2proc)
  - [CNVKIT](#cnvkit)
    - [pCNVkitAccess](#pcnvkitaccess)
    - [pCNVkitTarget](#pcnvkittarget)
    - [pCNVkitCov](#pcnvkitcov)
    - [pCNVkitRef](#pcnvkitref)
    - [pCNVkitFix](#pcnvkitfix)
    - [pCNVkitSeg](#pcnvkitseg)
    - [pCNVkitCall](#pcnvkitcall)
    - [pCNVkitPlot](#pcnvkitplot)
    - [pCNVkitRpt](#pcnvkitrpt)
    - [pCNVkit2Vcf](#pcnvkit2vcf)
  - [SEQ](#seq)
    - [pConsv](#pconsv)
    - [pGetPromoterBed](#pgetpromoterbed)
    - [pGetPromotersBed](#pgetpromotersbed)
  - [TCGA](#tcga)
    - [pDownload](#pdownload)
    - [pSample2SubmitterID](#psample2submitterid)
    - [pConvertExpFiles2Matrix](#pconvertexpfiles2matrix)
    - [pConvertMutFiles2Matrix](#pconvertmutfiles2matrix)
  - [WXSSTAT](#wxsstat)
    - [pVcf2List](#pvcf2list)
    - [pCallRate](#pcallrate-1)
    - [pCoverageByBamstats](#pcoveragebybamstats)
    - [pPlotBamstats](#pplotbamstats)
    - [pSnpEff2Stat](#psnpeff2stat)
    - [pPlotSnpEff](#pplotsnpeff)
  - [SNPARRAY](#snparray)
    - [pSNP6Genotype](#psnp6genotype)
    - [pGenoToAvInput](#pgenotoavinput)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->


# Documentation for bioprocs v0.0.1a
A set of procs for bioinformatics using [pyppl](https://github.com/pwwang/pyppl)

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
- - Both rownames and header of `full` and `triangle` can be omitted, just set `args.rownames = "NULL"` and `args.header = False`

#### output
- `outfile:file`: The output coordinate file

#### args
- `informat`: The format of the input file: full, triangle or pair. Default: full
- `rownames`: The `row.names` argument for `read.table`, default: 1
- `header`:   The `header` argument for `read.table` to read the input file, default: True.
- `k`:        How many dimension? Default: 2 (R^2)


###  pCluster
#### description
- Use `optCluster` to select the best number of clusters and cluster method, then perform the clustering

#### input
- `infile:file`: The input matrix file. Clustering will be performed against rows. If not, use `args.transpose` = True

#### output
- `outfile:file`: The output cluster file
- `outdir:dir`:   The output directory containing the figures

#### args
- `transpose`:    Transpose the input matrix. Default: False
- `header`:       Whether the input matrix contains header before transposing. Default: False
- `rownames`:     Which column is the rownames before transposing. Default: 1
- `plot`:         Whether plot the cluster. Default: True
- `nc`:           Number of clusters to test. Default: "2:15"
- `methods`:      The methods to test. Default: "all"
- - Could be any of "agnes", "clara", "diana", "fanny", "hierarchical", "kmeans", "model", "pam", "som", "sota", "em.nbinom", "da.nbinom", "sa.nbinom", "em.poisson", "da.poisson", "sa.poisson"
- - Multiple methods could be separated by comma (,), or put in a list
- - By default, fanny, model and sota will be excluded because fanny causes error and the latter two are slow. You can manually include them if you want.
- - Improper methods will be automatically excluded by `args.isCount`
- `isCount`:      Whether the data is count data. Corresponding methods will be tested. Default: False

#### requires
- [`r-optCluster`](https://rdrr.io/cran/optCluster/man/optCluster.html)
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
- `transpose`:Whether to transpose the matrix before cluster. Default: False

#### requires
- [`r-fastcluster`](https://cran.r-project.org/web/packages/fastcluster/index.html) if `args.fast` is True
- [`r-ggdendro`](https://cran.r-project.org/web/packages/ggdendro/index.html) if `args.gg` is True


## FASTX

###  pFastqPESim
#### description
- Simulate reads

#### input
- `in`: Index of the job/simulation, typically use range(10) for 10-time simulations

#### output
- `fq1:file`: The first pair read file
- `fq2:file`: The second pair read file

#### args
- `tool`:  The tool used for simulation. Default: wgsim (dwgsim)
- `len1`:  The length of first pair read. Default: 100
- `len2`:  The length of second pair read. Default: 100
- `num`:   The number of read PAIRs. Default: 1000000
- `seed`:  The seed for randomization. Default: None
- `gz`:    Whether generate gzipped read file. Default: True
- `wgsim`: The path of wgsim. Default: wgsim
- `dwgsim`:The path of wgsim. Default: dwgsim
- `ref`:   The reference genome. Required
- `params`:Other params for `tool`. Default: ""

#### requires
- [`wgsim`](https://github.com/lh3/wgsim)


###  pFastQC
#### description
- QC report for fastq file

#### input
- `fq:file`:    The fastq file (also fine with gzipped)

#### output
- `outdir:dir`: The output direcotry

#### args
- `tool`:    The tool used for simulation. Default: fastqc 
- `fastqc`:  The path of fastqc. Default: fastqc
- `nthread`: Number of threads to use. Default: 1
- `params`:Other params for `tool`. Default: ""

#### requires
- [`fastqc`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)


###  pFastMC
#### description
- Multi-QC based on pFastQC

#### input
- `qcdir:file`:  The direcotry containing QC files

#### output
- `outdir:dir`: The output direcotry

#### args
- `tool`:    The tool used for simulation. Default: multiqc 
- `multiqc`: The path of fastqc. Default: multiqc
- `params`:  Other params for `tool`. Default: ""

#### requires
- [`multiqc`](http://multiqc.info/)


###  pFastqPETrim
#### description
- Trim pair-end FASTQ reads

#### input
- `fq1:file`:  The input fastq file
- `fq2:file`:  The input fastq file

#### output
- `outfq1:file`: The trimmed fastq file
- `outfq2:file`: The trimmed fastq file

#### args
- `tool`        : The tools used for trimming. Default: trimmomatic (cutadapt|skewer)
- `cutadapt`    : The path of seqtk. Default: cutadapt
- `skewer`      : The path of fastx toolkit trimmer. Default: skewer
- `trimmomatic` : The path of trimmomatic. Default: trimmomatic
- `params`      : Other params for `tool`. Default: ""
- `nthread`     : Number of threads to be used. Default: 1
- - Not for cutadapt
- `gz`          : Whether gzip output files. Default: True
- `mem`         : The memory to be used. Default: 4G
- - Only for trimmomatic
- `minlen`      : Discard trimmed reads that are shorter than `minlen`. Default: 18
- - For trimmomatic, the number will be `minlen`*2 for MINLEN, as it filters before trimming
- `minq`        : Minimal mean qulity for 4-base window or leading/tailing reads. Default: 3
- `cut5`        : Remove the 5'end reads if they are below qulity. Default: 3
- `cut3`        : Remove the 3'end reads if they are below qulity. Default: 3
- - Not for skewer
- `adapter1`    : The adapter for sequence. Default: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
- `adapter2`    : The adapter for pair-end sequence. Default: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA

#### requires
- [`cutadapt`](http://cutadapt.readthedocs.io/en/stable/guide.html)
- [`skewer`](https://github.com/relipmoc/skewer)
- [`trimmomatic`](https://github.com/timflutre/trimmomatic)


###  pFastqSETrim
#### description
- Trim single-end FASTQ reads

#### input
- `fq:file`:  The input fastq file

#### output
- `outfq:file`: The trimmed fastq file

#### args
- `tool`        : The tools used for trimming. Default: trimmomatic (cutadapt|skewer)
- `cutadapt`    : The path of seqtk. Default: cutadapt
- `skewer`      : The path of fastx toolkit trimmer. Default: skewer
- `trimmomatic` : The path of trimmomatic. Default: trimmomatic
- `params`      : Other params for `tool`. Default: ""
- `nthread`     : Number of threads to be used. Default: 1
- - Not for cutadapt
- `gz`          : Whether gzip output files. Default: True
- `mem`         : The memory to be used. Default: 4G
- - Only for trimmomatic
- `minlen`      : Discard trimmed reads that are shorter than `minlen`. Default: 18
- - For trimmomatic, the number will be `minlen`*2 for MINLEN, as it filters before trimming
- `minq`        : Minimal mean qulity for 4-base window or leading/tailing reads. Default: 3
- `cut5`        : Remove the 5'end reads if they are below qulity. Default: 3
- `cut3`        : Remove the 3'end reads if they are below qulity. Default: 3
- - Not for skewer
- `adapter`     : The adapter for sequence. Default: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC

#### requires
- [`cutadapt`](http://cutadapt.readthedocs.io/en/stable/guide.html)
- [`skewer`](https://github.com/relipmoc/skewer)
- [`trimmomatic`](https://github.com/timflutre/trimmomatic)


###  pFastqSE2Sam
#### description
- Cleaned paired fastq (.fq, .fq.gz, .fastq, .fastq.gz file to mapped sam/bam file

#### args
- `tool`:   The tool used for alignment. Default: bwa (bowtie2|ngm)
- `bwa`:    Path of bwa, default: bwa
- `ngm`:    Path of ngm, default: ngm
- `bowtie2`:Path of bowtie2, default: bowtie2
- `rg`:     The read group. Default: {'id': '', 'pl': 'Illumina', 'pu': 'unit1', 'lb': 'lib1', 'sm': ''}
- - `id` will be parsed from filename with "_LX_" in it if not given
- - `sm` will be parsed from filename
- `ref`:    Path of reference file
- `params`: Other params for tool, default: ''


###  pFastqPE2Sam
#### description
- Cleaned paired fastq (.fq, .fq.gz, .fastq, .fastq.gz file to mapped sam/bam file

#### args
- `tool`   : The tool used for alignment. Default: bwa (bowtie2, ngm, star)
- `bwa`    : Path of bwa, default: bwa
- `ngm`    : Path of ngm, default: ngm
- `star`   : Path of ngm, default: STAR
- `bowtie2`: Path of bowtie2, default: bowtie2
- `rg`:     The read group. Default: {'id': '', 'pl': 'Illumina', 'pu': 'unit1', 'lb': 'lib1', 'sm': ''}
- - `id` will be parsed from filename with "_LX_" in it if not given
- - `sm` will be parsed from filename
- `ref`    : Path of reference file
- `refgene`: The GTF file for STAR to build index. It's not neccessary if index is already been built. Default: ''
- `params` : Other params for tool, default: ''


## WXSDOWN

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
- `outdir:dir`: The output directory

#### args
- `mutsig`: The path to `run_MutSigCV.sh`, default: 'mutsig'
- `mcr`: The Matlab MCR path

#### requires
- [MutSing](http://archive.broadinstitute.org/cancer/cga/mutsig_download)


###  pVcf2Maf
#### description
- Convert a snpEff-annotated somatic mutation vcf file (with normal and tumor samples) to [maf](https://wiki.nci.nih.gov/display/TCGA/Mutation+Annotation+Format+(MAF)+Specification) file

#### input
- `infile:file`: vcf file

#### output
- `outfile:file`: The maf file

#### args
- `vepdata`: The path of vep data. Default: "" (default data dir of vep)
- `vep`: The path of vep excutable. Default: "vep"
- `vcf2maf`: The path of vcf2maf excutable. Default: "vcf2maf.pl"
- `reffile`: The reference fasta file.
- `nthread`: The number of threads used by vep. Default: 1
- `filtervcf`: The filter vcf
- `params`: Other parameters for `vcf2maf.pl`, default: ""

#### requires
- [vcf2maf.py](https://github.com/mskcc/vcf2maf)


###  pMergeMafs
#### description
- Merge MAF files

#### input
- `indir:file`: The directory containing MAF files to be merged

#### output
- `outfile:file`: The merged MAF file


###  pMutsig4Plot
#### description
- Prepare somatic mutations for  plotting

#### input
- `msdir:file`:   The mutsig output directory

#### output
- `outfile:file`:  The file for plotting
```
	#PANEL: Somatic mutations
	#INFO: MT|PI
	#DESC: Mutation type|Putative impact
	# could also be bordercolor, you can have up to 4 shape features
	#TYPE: shape|bgcolor
	# could also be continuous
	# expressions for set: a,b,c
	#                 norminal: no
	#                 continuous: [0,1]
	#DATA: set|norminal
	#NCOL: 2|2
	#NAME_MT: Frameshift|Missense|Nonsense|Silent|Splice_site|TSS|Nonstop
	#NAME_PI: HIGH|MODERATE|LOW|MODIFIER
	#VALUE_MT: 0|1|20|13|4|17|14
	#EXP_MT: frameshift_variant,inframe_deletion,inframe_insertion|missense_variant,initiator_codon_variant,stop_retained_variant,rare_amino_acid_variant|stop_gained|synonymous_variant|splice_acceptor_variant,splice_donor_variant|start_lost,start_retained|stop_lost
	#
	Sample1	Sample2	Sample3	Sample4	Sample5
	ABC	missense_variant|HIGH	missense_variant|HIGH	...
	...
```

#### args
- `topn`:     the cutoff to select genes. If it is >= 1, top N genes will be selected, otherwise, it will be used as pvalue cutoff. Default: .05

#### requires
- [`pyvcf`](https://github.com/jamescasbon/PyVCF)


###  pMutPlot
#### description
- Plot mutations
```
	|           |             |           |           |---
	|- ftWidth -|  s   s   s  |- pnWidth -|- lgWidth -| snHeight
	|           |             |           |           |---
	    feature1
		feature2
```

#### input
- `indir:file`:    The input directory containing plot files

#### output
- `outfile:file`:  The plot png file


###  pCepip
#### input
- `avinput:file`: The avinput file
- `cell`:         The cell

#### output
- `outfile:file`: The cepip result file

#### args
- `bin-cepip`:    The jar file path of cepip, default: /data2/junwenwang/shared/tools/cepip/cepip.jar

#### requires
- [`cepip`](http://jjwanglab.org/cepip/)


## VCF

###  pVcfFilter
#### description
- Filter records in vcf file.

#### input
- `infile:file`: The input file

#### output
- `outfile:file`: The output file

#### args
- `tool`: Which tool to use for filtering. Default: 'vcflib'
- `vcflib_vcffilter`: The path of vcffilter from vcflib. Default: 'vcffilter'
- `gatk`            : The path of gatk. Default: 'gatk'
- `snpsift`         : The path of snpsift. Default: 'SnpSift'
- `bcftools`        : The path of bcftools. Default: 'bcftools'
- `samtools`        : The path of samtools. Default: 'samtools' (used by gatk to generate reference index)
- `picard`          : The path of picard. Default: 'picard' (used by picard to generate reference dict) 
- `params`          : Other params of `tool`. Default: ""
- `mem`             : The memory to be used. Default: "4G" (only for snpsift and gatk)
- `gz`              : Whether to gzip the output file. Default: False
- `keep`            : Whether to keep the filtered records. Default: True. (only for gatk, snpsift at filter step)
- `ref`             : The path of reference. Default: "" (for gatk)
- `tmpdir`          : The path of tmpdir. Default: <system tmpdir> (only used by gatk and snpsift)
- `nthread`         : The path of Default: 1	
- `selectors`:   Select records by:
- - type (snp, indel), sample genotypes (0, 1, 2), min genotype quality, filter (PASS, .)
- - for example:
```
		{"type": "snp", "genotype": {0: '0/0'}, "qual": 30}
		to select snps and whose genotype is '0/0' in 1st sample with quality >= 30
		{"genotype": {0: ['1/1', '0|1']}, "filter": ["PASS"]}
		to select records with PASS and genotype in 1st sample is '1/1' or '0/1'
```
- `filters`:     Filters depend on the tool you use on INFO filelds
- - format: `{"name1": "expression1", ...}`
- - If a string is specified, will convert to `{<tool name>: <expression>}`
- - Remember it filters OUT the records when ANY of the expression is true

#### requires
- [`pyvcf`](https://github.com/jamescasbon/PyVCF)
- [`gatk`](https://software.broadinstitute.org/gatk)
- [`bcftools`](http://www.htslib.org/doc/bcftools-1.2.html)
- [`snpsift`](http://snpeff.sourceforge.net/SnpSift.version_4_0.html)
- [`samtools`](https://github.com/samtools/samtools) if `args.ref` is not indexed, or bamutil is used for bam index file generation.
- [`picard`](https://broadinstitute.github.io/picard/command-line-overview.html) if `args.ref is not dicted.`


###  pVcfAnno
#### description
- Annotate the variants in vcf file.
- You have to prepare the databases for each tool.

#### input
- `infile:file`: The input vcf file

#### output
- `outfile:file`: The output file (output file of annovar will also be converted to vcf)
- `outdir`: The output directory, used to fetch some stat/summary files

#### args
- `tool`:            The tool used to do annotation. Default: snpeff
- `snpeff`:          The path of snpeff. Default: snpEff
- `vep`:             The path to vep. Default: vep
- `gz`:              Whether to gzip the result file. Default: False
- `annovar`:         The path of annovar. Default: annotate_variation.pl
- `annovar_convert`: The path of convert2annovar.pl, used to convert vcf to annovar input file. Default: convert2annovar.pl
- `genome`:          The genome for annotation. Default: hg19
- `tmpdir`:          The tmpdir, mainly used by snpeff. Default: <system tmpdir>
- `dbpath`:          The path of database for each tool. Required by 'annovar' and 'vep'
- `params`:          Other params for tool. Default: ''
- `snpeffStats`:     Whether to generate stats file when use snpeff. Default: False
- `mem`:             The memory used by snpeff. Default: '4G'

#### requires
- [`annovar`](http://doc-openbio.readthedocs.io/projects/annovar/en/latest/)
- [`snpeff`](http://snpeff.sourceforge.net/SnpEff_manual.html#intro)
- [`vep`](http://www.ensembl.org/info/docs/tools/vep/script/vep_tutorial.html)


###  pCallRate
#### description
- Calculate sample/snp call rate from single sample vcfs

#### input
- `indir:file`:     The dir containing the vcfs

#### output
- `outsample:file`: The report of call rate for each sample
- `figsample:file`: The bar chat of sample call rates
- `outsnp:file`:    The report of call rate for each snp
- `figsnp:file`:    The bar chat of snp call rates


## WXSCALL

###  pCNVnator
#### description
- Use `CNVnator` to call CNVs from bam file

#### input
- `infile:file`:  The bam file 

#### output
- `outfile:file`: The vcf file

#### args
- `cnvnator`:      The CNVnator executable, default: "cnvnator"
- `cnv2vcf`:  The converter executable to convert CNVnator results to vcf, default: "cnvnator2VCF.pl"
- `binsize`:  The bin_size, default: 100
- `genome`:   The genome: default: hg19
- `chrom`:    Chromosome names, default: "" (all chromosomes)
- `chrdir`:   The dir contains reference sequence of chromosomes, default: "" (don't specify)
- 

#### requires
- [CNVnator](https://github.com/abyzovlab/CNVnator)


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


###  pExpmat2Gct
#### description
- Convert expression matrix to GCT file.
- Refer to http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GCT for file format

#### input
- `expfile:file`: the input expression matrix file. Samples as columns, genes as rows.

#### output
- `outfile:file`: the gct file


###  pSampleinfo2Cls
#### description
- Convert sample infomation to cls file.
- Refer to http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#CLS for file format
- NOTE that the order of samples must be the same as in GMT file in further analysis.

#### input
- `sifile:file`: the sample information file.
- - Headers are: [Sample, ]Patient, Group, Batch
- - Rows are samples

#### output
- `outfile:file`: the cls file


###  pSSGSEA
#### description
- Single sample GSEA
- Refer to http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GCT for GCT file format
- Refer to http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GMT for GMT file format

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
- `nperm`:     Number of permutations. Default: 10000


###  pGSEA
#### description
- GSEA
- Refer to http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GCT for GCT file format
- Refer to http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#GMT for GMT file format
- Refer to http://software.broadinstitute.org/cancer/software/genepattern/file-formats-guide#CLS for CLS file format

#### input
- `gctfile:file`: the expression file
- `clsfile:file`: the class file
- `gmtfile:file`: the gmtfile for gene sets

#### output
- `outdir:file`: the output directory

#### args
- `weightexp`: Exponential weight employed in calculation of enrichment scores. Default: 0.75
- `nperm`:     Number of permutations. Default: 10000


###  pEnrichr
#### description
- Use APIs from http://amp.pharm.mssm.edu/Enrichr/help#api&q=1 to analyze a gene list

#### input
- `infile:file`: The gene list, each per line

#### output
- `outdir:dir`:  The output directory, containing the tables and figures.

#### args
- `topn`: Top N pathways used to plot. Default: 10
- `dbs`:  The databases to do enrichment against. Default: KEGG_2016
- - A full list can be found here: http://amp.pharm.mssm.edu/Enrichr/#stats
- - Multiple dbs separated by comma (,)
- `norm`: Normalize the gene list use [python-mygene](https://pypi.python.org/pypi/mygene/3.0.0)
- `rmtags`: Remove pathway tags in the plot. Default: True
- - For example: change "Lysine degradation_Homo sapiens_hsa00310" to "Lysine degradation".
- `plot`: Whether to plot the result. Default: True
- `title`: The title for the plot. Default: "Gene enrichment: {db}"

#### requires
- [python-mygene](https://pypi.python.org/pypi/mygene/3.0.0) if `args.norm` is `True`


###  pTargetEnrichr
#### description
- Use APIs from http://amp.pharm.mssm.edu/Enrichr/help#api&q=1 to analyze a gene list

#### input
- `infile:file`: The target genes with regulators
- - Format: 
- - Header is not required, but may specified in first line starting with `#`
- - If only 3 columns are there, the 3rd column is anyway the relation!
- - If only 4 columns are there, 3rd is target status, 4th is relation!
```
		  #Regulator	Target	Regulator status	Target status	Relation
		  has-mir-22	Gene	+	+	+
```

#### output
- `outdir:dir`:  The output directory, containing the tables and figures.

#### args
- `dbs`       : The databases to do enrichment against. Default: KEGG_2016
- - A full list can be found here: http://amp.pharm.mssm.edu/Enrichr/#stats
- - Multiple dbs separated by comma (,)
- `rmtags`    : Remove pathway tags in the plot. Default: True
- - For example: change "Lysine degradation_Homo sapiens_hsa00310" to "Lysine degradation".
- `enrplot`   : Whether to plot the result. Default: True
- `enrn`      : Top N pathways used to plot. Default: 10
- `netplot`   : Whether to plot the network. Default: True
- `netn`      : Top N pathways used to plot the network. Default: 5
- - Must <= `enrn`. If `netn` >= `enrn`, `netn` = `enrn`
- `title`     : The title for the plot. Default: "Gene enrichment: {db}"

#### requires
- [`python-mygene`](https://pypi.python.org/pypi/mygene/3.0.0) 
- [`graphviz`](https://pypi.python.org/pypi/graphviz)


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


## SNPINFO

###  pSNP2Bed
#### input
- `snpfile:file`: the snp file, each snp per line

#### output
- `outfile:file`: the result file, columns are:
- - chrom, start(0-based), end, name, score, strand, ref, allele

#### args
- `genome`: default: hg19
- `snpver`: default: snp147

#### requires
- [`python-cruzdb`](https://github.com/brentp/cruzdb)


###  pSNP2Avinput
#### input
- `snpfile:file`: the snp file, each snp per line

#### output
- `outfile:file`: the result avinput file

#### args
- `genome`: default: hg19
- `snpver`: default: snp147

#### requires
- [`python-cruzdb`](https://github.com/brentp/cruzdb)


## MARRAY

###  pCeldir2Matrix
#### description
- Convert CEL files to expression matrix
- File names will be used as sample names (colnames)

#### input
- `expdir:file`:  the directory containing the CEL files, could be gzipped

#### output
- `outfile:file`: the expression matrix file
- `outdir:dir`:   the directory containing expr file and plots

#### args
- `pattern` : The pattern to filter files. Default `'*'`
- `norm`    : The normalization method. Default: rma (mas5)
- `gfile`   : The group file. Default: ''
- `cdffile` : The cdffile. Default: ''
- `annofile`: The annotation file. Default: ''
- `exrows`  : Rows to be excluded, regular expression applied. Default: `[]`
- `boxplot` : Whether to plot a boxplot. Default: False
- `heatmap` : Whether to plot a heatmap. Default: False
- `histplot`: Whether to plot a histgram. Default: False
- `devpars` : Parameters for png. Default: `{'res': 300, 'width': 2000, 'height': 2000}`
- `boxplotggs`: The ggplot parameters for boxplot. Default: `['r:ylab("Expression")']`
- - See ggplot2 documentation.
- `heatmapggs`: The ggplot parameters for heatmap. Default: `['r:theme(axis.text.y = element_blank())']`
- `histplotggs`: The ggplot parameters for histgram. Default: `['r:labs(x = "Expression", y = "# Samples")']`	


###  aCelPat2Deg
#### description
- From celfils to degs with sample info file.

#### input
- `pattern`: The pattern to match the celfiles
- `sfile`  : The sample file


###  aCelPat2DegGSEA
#### description
- From celfils to degs with sample info file and do GSEA.

#### input
- `pattern`: The pattern to match the celfiles
- `sfile`  : The sample file
- `gmtkey` : The gmtkey to gmt file to do the GSEA. See `bioprocs.resource.pTxt`


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


## PLOT

###  pBoxplot
#### description
- Generate box plot

#### input
- `datafile:file`: The data file

#### output
- `outpng:file`: The output figure

#### args
- `header`:    `header` parameter for `read.table`, default: True
- `rownames`:  `row.names` parameter for `read.table`, default: 1
- `params`:    Other parameters for `boxplot`, default: ""


###  pScatterPlot
#### description
- Scatter plots with more information

#### input
- `infile:file`: The input file.
- - Format:
```
		X	Y	Size	Color
	A	1	1	1	1
	B	2	2	2	2
```
- - Column 3,4 can be omitted

#### output
- `outfile:file`: The plot

#### args
- `colfunc`: The functions to generate colors. Default: `heat.colors`
- - Available: rainbow, heat.colors, terrain.colors, topo.colors, and cm.colors.
- `type`:    The type of the symbols. Default: `circles`
- - Available: circles, squares, rectangles, stars, thermometers, boxplots.
- `inches`:  Scale the largest symbol to this size. Default: 1/3
- `data`:    The columns for render the symbols. Default: 0 (a simple dot plot)
- - circles:       3 (radii)
- - squares:       3 (length of sides)
- - rectangles:    3:4 (widths and heights)
- - starts:        3:? (?>5, a matrix with three or more columns giving the lengths of the rays from the center of the stars.)
- - thermometers:  3:? (?=5|6, The first two columns give the width and height of the thermometer symbols. If there are three columns, the third is taken as a proportion: the thermometers are filled (using colour fg) from their base to this proportion of their height. If there are four columns, the third and fourth columns are taken as proportions and the thermometers are filled between these two proportions of their heights. The part of the box not filled in fg will be filled in the background colour (default transparent) given by bg.) 
- - boxplots:      3:7 (a matrix with five columns. The first two columns give the width and height of the boxes, the next two columns give the lengths of the lower and upper whiskers and the fifth the proportion (with a warning if not in [0,1]) of the way up the box that the median line is drawn.)
- `main`:   The title of the plot. Default: NULL (the file name)
- `xlab`:   The labels for x axis. Default: colnames(mat)[1]
- `ylab`:   The labels for y axis. Default: colnames(mat)[2]
- `text`:   Whether to show the text of the symbols (rownames). Default: TRUE


###  pVenn
#### description
- Venn/UpsetR plots.

#### input
- `infile:file`: The input matrix
- - format:
```
			category1	category2	category3
		[e1]	0	1	1
		[e2]	0	0	1
		...
		[eN]	1	0	0
```
- rownames are not necessary but colnames are.

#### output
- `outfile:file`: The plot

#### args
- `tool`: Which tools to use. Default: auto (venn, upsetr, auto(n<=3: venn, otherwise upsetr))
- `rownames`: Whether input file has rownames. Default: False
- `vennParams`: Other params for `venn.diagram`. Default: {}
- `upsetParams`: Other params for `upset`. Default: {}

#### requires
- [`r-VennDiagram`](https://www.rdocumentation.org/packages/VennDiagram)
- [`r-UpSetR`](https://www.rdocumentation.org/packages/UpSetR)


## WXSANNO

###  pSnpEff
#### description
- This is the default command. It is used for annotating variant filed (e.g. VCF files).

#### input
- `infile:file`:  The input file 

#### output
- `outdir:file`: The directory containing output anntated file, snpEff_genes.txt and snpEff_summary.html

#### args
- `snpEff`:       The snpEff executable, default: "snpEff"
- `params`:    Other parameters for `snpEff`, default: "-Xms1g -Xmx4g -v"
- `genome`:    The genome used for annotation, default: "hg19"
- `informat`:  The format of input file [vcf or bed], default: "vcf"
- `outformat`: The format of output file [vcf, gatk, bed, bedAnn], default: "vcf"
- `csvStats`:  Whether to generate csv stats file, default: True.
- `htmlStats`: Whether to generate the html summary file, default: False.
- `javamem`:   The memory to use. Default: '-Xms1g -Xmx8g'

#### requires
- [snpEff](http://snpeff.sourceforge.net/SnpEff_manual.html)


## RNASEQ

###  pExpdir2Matrix
#### description
- Convert expression files to expression matrix
- File names will be used as sample names (colnames)
- Each gene and its expression per line.
- Suppose each expression file has the same rownames and in the same order.

#### input
- `expdir:file`:  the directory containing the expression files, could be gzipped

#### output
- `outfile:file`: the expression matrix file
- `outdir:dir`:   the directory containing expr file and plots

#### args
- `pattern` : The pattern to filter files. Default `'*'`
- `header`  : Whether each expression file contains header. Default: `False`
- `exrows`  : Rows to be excluded, regular expression applied. Default: `["^Sample", "^Composite", "^__"]`
- `boxplot` : Whether to plot a boxplot. Default: False
- `heatmap` : Whether to plot a heatmap. Default: False
- `histplot`: Whether to plot a histgram. Default: False
- `devpars` : Parameters for png. Default: `{'res': 300, 'width': 2000, 'height': 2000}`
- `boxplotggs`: The ggplot parameters for boxplot. Default: `['r:ylab("Expression")']`
- - See ggplot2 documentation.
- `heatmapggs`: The ggplot parameters for heatmap. Default: `['r:theme(axis.text.y = element_blank())']`
- `histplotggs`: The ggplot parameters for histgram. Default: `['r:labs(x = "Expression", y = "# Samples")']`	


###  pBatchEffect
#### description
- Remove batch effect with sva-combat.

#### input
- `expr:file`:  The expression file, generated by pExpdir2Matrix
- `batch:file`: The batch file defines samples and batches.

#### output
- `outfile:file`: the expression matrix file
- `outdir:dir`:   the directory containing expr file and plots

#### args
- `tool`    : The tool used to remove batch effect. Default `'combat'`
- `boxplot` : Whether to plot a boxplot. Default: False
- `heatmap` : Whether to plot a heatmap. Default: False
- `histplot`: Whether to plot a histgram. Default: False
- `devpars` : Parameters for png. Default: `{'res': 300, 'width': 2000, 'height': 2000}`
- `boxplotggs`: The ggplot parameters for boxplot. Default: `['r:ylab("Expression")']`
- - See ggplot2 documentation.
- `heatmapggs`: The ggplot parameters for heatmap. Default: `['r:theme(axis.text.y = element_blank())']`
- `histplotggs`: The ggplot parameters for histgram. Default: `['r:labs(x = "Expression", y = "# Samples")']`	


###  pRawCounts2
#### description
- Convert raw counts to another unit

#### input
- `expfile:file`: the expression matrix
- - rows are genes, columns are samples

#### output
- `outfile:file`: the converted expression matrix

#### args
- `transpose`: transpose the input matrix? default: False
- `log2`:      whether to take log2? default: False
- `unit`:      convert to which unit? default: cpm (or rpkm, tmm)
- `header`:    whether input file has header? default: True
- `rownames`:  the index of the column as rownames. default: 1
- `glenfile`:  the gene length file, for RPKM
- - no head, row names are genes, have to be exact the same order and length as the rownames of expfile
- `boxplot` : Whether to plot a boxplot. Default: False
- `heatmap` : Whether to plot a heatmap. Default: False
- `histplot`: Whether to plot a histgram. Default: False
- `devpars` : Parameters for png. Default: `{'res': 300, 'width': 2000, 'height': 2000}`
- `boxplotggs`: The ggplot parameters for boxplot. Default: `['r:ylab("Expression")']`
- - See ggplot2 documentation.
- `heatmapggs`: The ggplot parameters for heatmap. Default: `['r:theme(axis.text.y = element_blank())']`
- `histplotggs`: The ggplot parameters for histgram. Default: `['r:labs(x = "Expression", y = "# Samples")']`	

#### requires
- [edgeR](https://bioconductor.org/packages/release/bioc/html/edger.html) if cpm or rpkm is chosen
- [coseq](https://rdrr.io/rforge/coseq/man/transform_RNAseq.html) if tmm is chosen


###  pDeg
#### description
- Detect DEGs for RNA-seq data

#### input
- `efile:file`: The expression matrix
- `gfile:file`: The group information
- - Like:
```
		Sample1	Group1
		Sample2	Group1
		Sample3	Group1
		Sample4	group2
		Sample5	group2
		Sample6	group2
```

#### output
- `outfile:file`: The DEG list
- `outdir:file`:  The output directory containing deg list and plots

#### args
- `tool`      : the tool used to detect DEGs. Default: 'edger' (deseq2)
- `filter`    : filter out low count records. Default: `"1,2"` (At least 2 samples have at least 2 reads)
- `mdsplot`   : whether to plot the MDS plot, default : True
- `volplot`   : whether to plot the volcano plot, default : True
- `maplot`    : whether to plot MA plots within each group, default : False
- `heatmap`   : whether to plot the heatmap using DEGs. Default : False
- `heatmapn`  : How many genes to be used for heatmap. If `heatmapn`, the number will be `heatmapn * # DEGs`. Default: 100
- `heatmapggs`: The ggplots options for heatmap. Default : []
- `maplotggs` : The ggplots options for maplot. Default : []
- `volplotggs`: The ggplots options for volplot. Default : []
- `devpars`   : Parameters for png. Default: `{'res': 300, 'width': 2000, 'height': 2000}`


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

#### output
- `outfile:file`: The output file


## BEDTOOLS

###  pGetfasta
#### description
- `bedtools getfasta` extracts sequences from a FASTA file for each of the intervals defined in a BED/GFF/VCF file.

#### input
- `infile:file`: The input bed file
- `fafile:file`: The input fasta file

#### brings
- `fafile`: "{{fafile | fn}}.fa*i", The fasta index file

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


## RANK

###  pRankProduct
#### description
- Calculate the rank product of a set of ranks. Refer to [here](https://en.wikipedia.org/wiki/Rank_product)

#### input
- `infile:file`: The input file
- - Format:
```
				Case1	Case2	...
	Feature1	8.2  	10.1 	...
	Feature2	2.3  	8.0  	...
	...
```
- - Or instead of values, you can also have ranks in the input file:
```
				Rank1	Rank2	...
	Feature1	2    	1    	...
	Feature2	3    	2    	...
	...
```

#### output
- `outfile:file`: The output file with original ranks, rank products and p-value if required

#### args
- `informat`: The input format of the values. Whether they are real values (value) or ranks (rank). Default: value
- `pval`:     Whether to calculate the p-value or not. Default: True
- `header`:   Whether the input file has headers (rownames are required!). Default: True
- `plot`:     Number of rows to plot. Default: 0 (Don't plot)
- `cex`:      Font size for plotting. Default: 0.9
- `cnheight`: Colname height. Default: 80
- `rnwidth`:  Rowname width. Default: 50
- `width`:    Width of the png file. Default: 2000
- `height`:   height of the png file. Default: 2000


## WXSPREP

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
- `trimmomatic`:    The trimmomatic executable, default: "trimmomatic"
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
- `trimmomatic`:    The trimmomatic executable, default: "trimmomatic"
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
- `reffile:file`: The reference file

#### output
- `outfile:file`: The output sam file

#### args
- `bwa`:    The bwa executable, default: bwa
- `params`: Other params for bwa mem, default: "-M"
- `nthread`: 1

#### requires
- [bwa](https://github.com/lh3/bwa)


###  pAlignSEByBWA
#### description
- Align paired-end reads to reference genome using bwa mem

#### input
- `infile:file`:  read file (fastq, or fastq gzipped)
- `reffile:file`: The reference file

#### brings
- `reffile#bwt`: "{{reffile | bn}}.bwt", 
- `reffile#sa`:  "{{reffile | bn}}.sa",
- `reffile#ann`: "{{reffile | bn}}.ann",
- `reffile#amb`: "{{reffile | bn}}.amb",
- `reffile#pac`: "{{reffile | bn}}.pac"

#### output
- `outfile:file`: The output sam file

#### args
- `bwa`:    The bwa executable, default: bwa
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
- `reffile:file`: The reference file

#### output
- `outfile:file`: The output sam/bam file

#### args
- `ngm`:    The NextGenMap executable, default: ngm
- `nthread`: 1
- `outtype`: sam or bam, default: sam (only sam for now, due to bug of ngm 0.5.3 (fixed in 0.5.4))
- `params`: Other params for ngm, default: "--rg-id ngm --rg-sm sample"

#### requires
- [NextGenMap](https://github.com/Cibiv/NextGenMap/wiki)


###  pAlignSEByNGM
#### description
- Align single-end reads to reference genome using NextGenMap

#### input
- `infile1:file`: read file 1 (fastq, or fastq gzipped)
- `infile2:file`: read file 2 (fastq, or fastq gzipped)
- `reffile:file`: The reference file

#### output
- `outfile:file`: The output sam/bam file

#### args
- `ngm`:    The NextGenMap executable, default: ngm
- `nthread`: 1
- `outtype`: sam or bam, default: sam (only sam for now, due to bug of ngm 0.5.3 (fixed in 0.5.4))
- `params`: Other params for ngm, default: "--rg-id ngm --rg-sm sample"

#### requires
- [NextGenMap](https://github.com/Cibiv/NextGenMap/wiki)


###  pMergeBams
#### description
- Merge bam files

#### input
- `bamdir:dir`:   the dir containing bam files 

#### output
- `outfile:file`: the merged bam file

#### args
- `samtools`: the executable path of samtools, default: "samtools"
- `nthread`:      Number of BAM/CRAM compression threads
- `params`:       Other parameters for `samtools merge`, default: ""

#### requires
- [samtools](http://www.htslib.org/)


## DEG

###  pExpdirMatrix
#### description
- Convert expression files to expression matrix
- File names will be used as sample names (colnames)
- Each gene and its expression per line.

#### input
- `expdir:file`:  the directory containing the expression files, could be gzipped

#### output
- `expfile:file`: the expression matrix


###  pDEGByEdgeR
#### description
- Call DEG from expression matrix

#### input
- `expfile:file`: the expression matrix
- `group1`:       columns of group1 (separated by comma)
- `group2`:       columns of group2 (separated by comma)
- `group1name`:   the name of group1
- `group2name`:   the name of group2   

#### output
- `degdir:dir`:   the output directory containing DEGs and plots

#### args
- `filter`:  the pair (X,Y) on how to filter the data by cpm (`d <- d[rowSums(cpm(d)>X) >= Y,]`). Default: "1,2"
- - keep genes with at least X counts per million (CPM) in at least Y samples
- `pval`:    the cutoff of DEGs (default: .05)
- `paired`:  whether the samples are paired, default: False
- `bcvplot`: whether to plot biological coefficient of variation, default: True
- `displot`: whether to plot biological coefficient of variation, default: True
- `fcplot`:  whether to plot fold changes, default: True

#### requires
- [edgeR](https://bioconductor.org/packages/release/bioc/html/edger.html)


###  pMArrayLimma
#### description
- Call degs of microarray data by limma

#### input
- `expfile:file`: The expression matrix
- `group1`:      The 1st group
- `group2`:      The 2nd group
- `group1name`:  The name of 1st group
- `group2name`:  The name of 2nd group

#### output
- `degdir:dir`:   the output directory containing DEGs and plots

#### args
- `norm`:    the normalization methods, separated by comma. Support normalization methods: `quan` (quantile). Default: "quan"
- `boxplot`: draw boxplot? Default: True
- `paired`:  whether the samples are paired, default: False
- `filter`:  the pair (X,Y) on how to filter the data by expression (`d <- d[rowSums(d>X) >= Y,]`). Default: "1,2"
- - keep genes with at least X exp in at least Y samples
- `pval`:    the pvalue cutoff of DEGs (default: .05)
- `qval`:    the qvalue cutoff of DEGs (default: .05)
- `heatmap`: whether to plot heatmap, default: True
- `hmn`:     Number of gene used for heatmap, default: 50
- `hmmar`:   Margins for heatmap, default: "10,7"
- `volplot`: whether to plot the volcano plot, default: True

#### requires
- [limma](https://bioconductor.org/packages/release/bioc/html/limma.html)
- [ggplot2](https://bioconductor.org/packages/release/bioc/html/ggplot2.html)


###  pRawCounts2
#### description
- Convert raw counts to another unit

#### input
- `expfile:file`: the expression matrix
- - rows are genes, columns are samples, if not use `args.transpose = True`

#### output
- `outfile:file`: the converted expression matrix

#### args
- `transpose`: transpose the input matrix? default: False
- `log2`:      whether to take log2? default: False
- `unit`:      convert to which unit? default: cpm (or rpkm, tmm)
- `header`:    whether input file has header? default: True
- `rownames`:  the index of the column as rownames. default: 1
- `glenfile`:  the gene length file, for RPKM
- - no head, row names are genes, have to be exact the same order and length as the rownames of expfile

#### requires
- [edgeR](https://bioconductor.org/packages/release/bioc/html/edger.html) if cpm or rpkm is chosen
- [coseq](https://rdrr.io/rforge/coseq/man/transform_RNAseq.html) if tmm is chosen


## VCFNEXT

###  pStats2Matrix
#### description
- Convert csvstat file from snpEff to R-readable matrix for plotting

#### input
- `indir:file`: The directory containing the csv stat files from `snpEff ann`

#### output
- `outdir:dir`: The output directory

#### args
- `chroms`:     The chromsome filter. Default: "" (all chroms)
- - Note: snpEff csvstat file has no "chr" prefix


###  pPlotStats
#### description
- Plot snpEff annotation statistics

#### input
- `indir:file`: The snpEff result directory containing matrix files generated by pSnpEff2Stat

#### output
- `outdir:dir`: The output directory

#### requires
- [`pwwang/corrplot`](https://github.com/pwwang/corrplot)
- - use `library(devtools); install.github("pwwang/corrplot")`
- [`ggplot2`](http://ggplot2.org/)


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
- `picard`:     The picard executable, default: "picard"
- `params`:  Other parameters for picard MarkDuplicates, default: ""
- `tmpdir`:  The tmpdir to use. Default: /tmp

#### requires
- [picard](https://broadinstitute.github.io/picard/)


###  pAddOrReplaceReadGroups
#### description
- Replace read groups in a BAM file.This tool enables the user to replace all read groups in the INPUT file with a single new read group and assign all reads to this read group in the OUTPUT BAM file.
- For more information about read groups, see the [GATK Dictionary entry](https://www.broadinstitute.org/gatk/guide/article?id=6472). 
- This tool accepts INPUT BAM and SAM files or URLs from the Global Alliance for Genomics and Health (GA4GH) (see http://ga4gh.org/#/documentation).

#### input
- `infile:file`:  The bam file
- `rg`:           The read group information. For example:
- - "RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM=20"

#### output
- `outfile:file`: The bam file with read group added

#### args
- `picard`:     The picard executable, default: "picard "
- `params`:  Other parameters for picard AddOrReplaceReadGroups, default: ""

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
- `picard`:     The picard executable, default: "picard"
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
- `picard`:     The picard executable, default: "picard"
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
- `picard`:     The picard executable, default: "picard"
- `order`:   The sort order, default: coordinate. Possible: unsorted, queryname, coordinate, duplicate
- `outtype`: The type of output file, sam or bam. Default: bam
- `params`:  Other parameters for `picard SortSam`, default: ""
- `tmpdir`:  The tmpdir to use. Default: /tmp
- `javamem`: The memory for java vm. Default: "-Xms1g -Xmx8g"

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
- `picard`:    The picard executable, default: "picard"
- `params`:  Other parameters for `picard BuildBamIndex`, default: "-Xms1g -Xmx8g"

#### requires
- [picard](http://broadinstitute.github.io/picard/command-line-overview.html)


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


## BED

###  pBedSort
#### description
- Sort bed files

#### input
- `infile:file`: The input file

#### output
- `outfile:file`: The output file

#### args
- `tool`:         The tool used to sort the file. Default: sort (bedtools, bedops)
- `bedtools`:     The path to bedtools. Default: bedtools
- `bedops_sort`:  The path to bedops' sort-bed. Default: sort-bed
- `sort`:         The path to linux's sort. Default: sort
- `mem`:          The memory to use. Default: 8G
- `tmpdir`:       The tmpdir to use. Default: `$TMPDIR`
- `unique`:       Remove the dupliated records? Default: True
- `params`:       Other params for `tool`. Default: ''

#### requires
- [`bedtools`](http://bedtools.readthedocs.io/en/latest/index.html)
- [`bedops`](https://github.com/bedops/bedops)


###  pBedIntersect
#### description
- Find intersections of two bed files.
- Input files must be sorted.

#### input
- `infile1:file`: The 1st input bed file
- `infile2:file`: The 2nd input bed file

#### output
- `outfile:file`: The output file

#### args
- `tool`:         The tool used to sort the file. Default: bedtools (bedops)
- `bedtools`:     The path to bedtools. Default: bedtools
- `bedops`:  The path to bedops. Default: bedops
- `params`:       Other params for `tool`. Default: ''

#### requires
- [`bedtools`](http://bedtools.readthedocs.io/en/latest/index.html)
- [`bedops`](https://github.com/bedops/bedops)


###  pBedCluster
#### description
- Assign cluster id to each record

#### input
- `infile:file`: The input bed file

#### output
- `outfile:file`: The output file

#### args
- `tool`:         The tool used to sort the file. Default: bedtools
- `bedtools`:     The path to bedtools. Default: bedtools
- `params`:       Other params for `tool`. Default: ''

#### requires
- [`bedtools`](http://bedtools.readthedocs.io/en/latest/index.html)


## RESOURCE

###  pTxt
#### description
- Download CSV format files.

#### input
- `in`: The name of the resource

#### output
- `outfile:file`: The output file

#### args
- `cols`:      Select the columns to keep. Default: '' (all cols)
- `rowfilter`: Filter rows. For example, to filter out rows not start with 'Chr':
- - `"lambda x: not x[0].startswith('Chr')"`
- - Note that rowfilter applied before cols filter.
- `urls`:      Available resources and their urls.
- `gz`:        Whether to gzip the output file.

#### requires
- [`curl`](https://en.wikipedia.org/wiki/CURL)


## PCA

###  pPCA
#### description
- Perform PCA analysis

#### input
- `infile:file`: The matrix to do the analysis
- - Note that rows are samples, columns are features, if not, use `args.transpose = True`

#### output
- `outfile:file`: The output coordinate file
- - Columns are PCs, rows are samples

#### args
- `transpose`:  Whether to transpose the input matrix from infile. Default: False
- `rownames`:   The `row.names` argument for `read.table`, default: 1
- `header`:     The `header` argument for `read.table` to read the input file, default: True.
- `screeplot`:  Whether to generate the screeplot or not. Default: True
- `sp_ncp`:     Number of components in screeplot. Default: 0 (auto detect)
- - if total # components (tcp) < 20: use all
- - else if tcp > 20, use 20
- `varplot`:    Whether to generate the variable plot or not. Default: False
- `biplot`:     Whether to generate the variable plot or not. Default: True

#### requires
- [`r-factoextra`](https://cran.r-project.org/web/packages/factoextra/index.html) for plots


###  pSelectPCs
#### description
- Select a subset of PCs from pPCA results

#### input
- `indir:file`: The directory generated from pPCA

#### output
- `outfile:file`: The file containing selected PCs

#### args
- `n`: The number of PCs to select. Default: 0.9
- - If it is < 1, used as the % variation explained from stdev.txt


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


## SAMBAM

###  pSam2Bam
#### description
- Deal with mapped sam/bam files, including sort, markdup, and/or index

#### input
- `infile:file`: The input file

#### output
- `outfile:file`: The output bam file
- `idxfile:file`: The index of the output bam file
- - If args.index == False, it'll a link to outfile and should be never used

#### args
- `tool`             : The tool used to do the sort. Default: sambamba (picard|sambamba|biobambam|samtools)
- `sambamba`         : The path of the sambamba. Default: sambamba 
- `picard`           : The path of the picard. Default: picard 
- `biobambam_bamsort`: The path of the biobambam's bamsort. Default: bamsort 
- `samtools`         : The path of the samtools. Default: samtools 
- `sort`             : Do sorting? Default: True 
- - If input is sam, tool is biobambam, this should be True
- `index`            : Do indexing? Default: True
- `markdup`          : Do duplicates marking? Default: False
- - `rmdup` for samtools will be called
- `rmdup`            : Do duplicates removing? Default: False
- `tmpdir`           : The tmp dir used to store tmp files. Default: <system default tmpdir>
- `sortby`           : Sort by coordinate or queryname. Default: coordinate
- `nthread`          : Default: 1
- `informat`         : The format of input file. Default: <detect from extension> (sam|bam)
- `params`           : Other parameters for `tool`. Defaut: ""
- `mem`              : The max memory to use. Default: "16G"
- - Unit could be G/g/M/m
- - Will be converted to -Xmx4G, and -Xms will be 1/8 of it

#### requires
- [sambamba](https://lomereiter.github.io/sambamba/docs/sambamba-view.html) if `args.tool` == samtools or reference used but not indexed.
- [picard](https://broadinstitute.github.io/picard/command-line-overview.html)
- [biobambam](https://github.com/gt1/biobambam2)
- [samtools](https://github.com/samtools/samtools)


###  pBamMarkdup
#### description
- Mark/remove duplicates for bam files

#### input
- `infile:file`: The input file

#### output
- `outfile:file`: The output bam file

#### args
- `tool`             : The tool used to do the sort. Default: sambamba (picard|sambamba|biobambam|samtools|bamutil)
- `sambamba`         : The path of sambamba. Default: sambamba 
- `picard`           : The path of picard. Default: picard 
- `biobambam_bamsort`: The path of biobambam's bamsort. Default: bamsort 
- `samtools`         : The path of samtools. Default: samtools 
- `bamutil`          : The path of bamutil. Default: bam
- `rmdup`            : Do duplicates removing? Default: False
- - Samtools will anyway remove the duplicates
- `tmpdir`           : The tmp dir used to store tmp files. Default: <system default tmpdir>
- `nthread`          : Default: 1
- - Not available for samtools and picard
- `params`           : Other parameters for `tool`. Defaut: ""
- `mem`              : The max memory to use. Default: "16G"
- - Unit could be G/g/M/m
- - Will be converted to -Xmx4G, and -Xms will be 1/8 of it

#### requires
- [sambamba](https://lomereiter.github.io/sambamba/docs/sambamba-view.html)
- [picard](https://broadinstitute.github.io/picard/command-line-overview.html)
- [biobambam](https://github.com/gt1/biobambam2)
- [samtools](https://github.com/samtools/samtools)
- [bamutil](http://genome.sph.umich.edu/wiki/BamUtil#Programs)


###  pBamRecal
#### description
- Recalibrate a bam file

#### input
- `infile:file`: The bam file

#### brings
- `infile`: {{in.infile | bn}}.bai, the index file of bam

#### output
- `outfile:file`: The output bam file

#### args
- `tool`                         : The tool used to recalibrate the bam file. Default: `gatk` (gatk|bamutil)
- `gatk`                         : The path of gatk, including java path. Default: `gatk`
- `samtools`                     : The path of samtools. Default: `samtools`
- `bamutil`                      : The path of bamutil. Default: `bam`
- `picard`                       : The path of picard. Default: `picard`
- `paramsRealignerTargetCreator` : Other parameters for `gatk RealignerTargetCreator`. Defaut: ""
- `paramsIndelRealigner`         : Other parameters for `gatk IndelRealigner`. Defaut: ""
- `paramsBaseRecalibrator`       : Other parameters for `gatk BaseRecalibrator`. Defaut: ""
- `paramsPrintReads`             : Other parameters for `gatk PrintReads`. Defaut: ""
- `params`                       : Other parameters for `bam recab`. Default: ""
- `mem`                          : The max memory to use. Default: "32G"
- `knownSites`                   : The known polymorphic sites to mask out. Default: "" (Required for GATK)
- `ref`                          : The reference file. Required.
- - Will be converted to -Xmx4G, and -Xms will be 1/8 of it

#### requires
- [gatk](https://software.broadinstitute.org/gatk)
- [samtools](https://github.com/samtools/samtools) if `args.ref` is not indexed, or bamutil is used for bam index file generation.
- [picard](https://broadinstitute.github.io/picard/command-line-overview.html) if `args.ref is not dicted.`


###  pBamReadGroup
#### description
- Add or replace read groups of a bam file

#### input
- `infile:file`: The bam file

#### output
- `outfile:file`: The output bam file

#### args
- `tool`                         : The tool used. Default: `picard` (picard|bamutil)
- `picard`                       : The path of picard. Default: `picard`
- `bamutil`                      : The path of bamutil. Default: `bam`
- `rg`                           : The read group. Default: {'id': '', 'pl': 'Illumina', 'pu': 'unit1', 'lb': 'lib1', 'sm': ''}
- - `id` will be parsed from filename with "_LX_" in it if not given
- - `sm` will be parsed from filename
- `params`                       : Other parameters for `tool`. Defaut: ""
- `mem`                          : The max memory to use. Default: "4G"
- - Will be converted to -Xmx4G, and -Xms will be 1/8 of it
- `tmpdir`                       : The temporary directory. Default: <system tmpdir>

#### requires
- [gatk](https://lomereiter.github.io/sambamba/docs/sambamba-view.html)
- [samtools](https://github.com/samtools/samtools) if `args.ref` is not indexed.
- [picard](https://broadinstitute.github.io/picard/command-line-overview.html) if `args.ref is not dicted.`


###  pBamReorder
#### description
- Reorder a sam/bam file by a given reference file using `picard ReorderSam`

#### input
- `infile:file`: The sam/bam file

#### output
- `outfile:file`: The output bam file

#### args
- `picard`                       : The path of picard. Default: `picard`
- `ref`                          : The reference file. Required
- `params`                       : Other parameters for `picard ReorderSam`. Defaut: ""
- `mem`                          : The max memory to use. Default: "4G"
- - Will be converted to -Xmx4G, and -Xms will be 1/8 of it
- `tmpdir`                       : The temporary directory. Default: <system tmpdir>

#### requires
- [picard](https://broadinstitute.github.io/picard/command-line-overview.html)


###  pBamMerge
#### description
- Merges multiple SAM and/or BAM files (must be sorted by coordinate) into a single file.

#### input
- `inlist:file`: The directory containing sam/bam files to be merged

#### output
- `outfile:file`: The merged bam file

#### args
- `tool`     : The tool used to merge. Default: bamutil (picard|samtools|sambamba)
- `picard`   : The path of picard. Default: `picard`
- `bamutil`  : The path of bamutil. Default: `bam`
- `samtools` : The path of samtools. Default: `samtools`
- `sambamba` : The path of sambamba. Default: `sambamba`
- `params`   : Other parameters for `tool`. Defaut: ""
- `mem`      : The max memory to use. Default: "4G"
- - Will be converted to -Xmx4G, and -Xms will be 1/8 of it, just for picard
- `tmpdir`   : The temporary directory. Default: <system tmpdir>
- `nthread`  : # threads to use. Default: 1
- - For picard, if nthread>1, USE_THREADING=true, otherwise USE_THREADING=false

#### requires
- [picard](https://broadinstitute.github.io/picard/command-line-overview.html)


###  pBam2Gmut
#### description
- Call germline (snps and indels) from a call-ready bam file.

#### input
- `infile:file`: The input bam file

#### brings
- `infile`: `{{in.infile | bn}}.bai`, the bam index file

#### output
- `outfile:file`: The vcf file containing the mutations

#### args
- `tool`:         The tool used to call mutations. Default: gatk (vardict, snvsniffer, platypus, strelka)
- `gatk`:         The path of gatk. Default: gatk
- `vardict`:      The path of vardict. Default: vardict
- `snvsniffer`:   The path of snvsniffer. Default: SNVSniffer
- `samtools`:     The path of samtools. Default: samtools (used to generate reference index)
- `platypus`:     The path of platypus. Default: platypus
- `strelka`:      The path of strelka. Default: configureStrelkaGermlineWorkflow.py
- `configParams`: The params for `strelka` configuration. Default: ""
- `picard`:       The path of picard. Default: picard
- `mem`:          The memory to be used. Default: 32G
- - will be converted to -Xms4G -Xmx32G for java programs
- `ref`:          The reference file. Required.
- `gz`:           Gzip output file? Default: False
- `tmpdir`:       The temporary directory. Default: <system tmpdir>
- `params`:       Other params for `tool`. Default: ""

#### requires
- [gatk](https://lomereiter.github.io/sambamba/docs/sambamba-view.html)
- [samtools](https://github.com/samtools/samtools) if `args.ref` is not indexed.
- [picard](https://broadinstitute.github.io/picard/command-line-overview.html) if `args.ref is not dicted.`
- [vardict](https://github.com/AstraZeneca-NGS/VarDict)
- [snvsniffer](http://snvsniffer.sourceforge.net/homepage.htm#latest)
- [platypus](http://www.well.ox.ac.uk/platypus)
- [strelka@2.7.1+](https://github.com/Illumina/strelka)


###  pBam2Cnv
#### description
- Detect copy number variation from bam files.

#### input
- `input:file`: The bam file

#### brings
- `infile`: "{{in.infile | bn}}.bai" The bam index file

#### output
- `outfile:file`: The output vcf file
- `outdir`: The output directory containing other result files

#### args
- `gz`                    : Whether to gzip the output vcf file. Default: False
- `tool`                  : The tool used to call cnv. Default: 'cnvkit'
- `cnvnator`              : The path of cnvnator. Default: 'cnvnator'
- `cnvnator2vcf`          : The path of cnvnator2VCF. Default: 'cnvnator2VCF.pl'
- `cnvkit`                : The path of cnvkit. Default: 'cnvkit.py'
- `wandy`                 : Tha path of Wandy. Default: 'Wandy'. A `tool.info` file should be with the executable file.
- `ref`                   : The reference file. Required by cnvkit to generate access file. Default: ''
- `cnvkitAccessParams`    : The params for cnvkit access command. Default: '-s 5000'
- `cnvkitTargetParams`    : The params for cnvkit target command. Default: '--split --short-names'
- `cnvkitCoverageParams`  : The params for cnvkit coverage command. Default: ''
- `cnvkitReferenceParams` : The params for cnvkit reference command. Default: '--no-edge'
- `cnvkitFixParams`       : The params for cnvkit fix command. Default: '--no-edge'
- `cnvkitSegmentParams`   : The params for cnvkit segment command. Default: ''
- `cnvkitCallParams`      : The params for cnvkit call command. Default: ''
- `cnvkitPlotParams`      : The params for cnvkit plot command. Default: ''
- `cnvkitBreaksParams`    : The params for cnvkit breaks command. Default: ''
- `cnvkitGainlossParams`  : The params for cnvkit gainloss command. Default: ''
- `cnvkitMetricsParams`   : The params for cnvkit metrics command. Default: ''
- `cnvkitSegmetricsParams`: The params for cnvkit segmetrics command. Default: '--iqr'
- `cnvkitExportParams`    : The params for cnvkit export command. Default: ''
- `cnvkitScatterParams`   : The params for cnvkit scatter command. Default: [''] # multiple scatter plots
- `cnvkitHeatmapParams`   : The params for cnvkit heatmap command. Default: [''] # multiple heatmap plots
- `cnvkitDiagramParams`   : The params for cnvkit diagram command. Default: ''
- `cnvkitReport`          : Generate cnvkit reports? Default: True
- `cnvkitPlot`            : Generate cnvkit plots? Default: True
- `cnvnatorBinsize`       : Bin size for cnvnator. Default: 100
- `cnvnatorGenome`        : Genome for cnvnator. Default: 'hg19'. (NCBI36, hg18, GRCh37, hg19)
- `params`                : The params for `tool`. Default: '-t 1' # wandy 1:hg19 solid cell/blood, 2:hg19 cell free/plamsa, 3:hg38 solid cell/blood, 4:hg38 cell free/plamsa
- `mem`                   : The memory used. Default: '20G' # only for wandy
- `nthread`               : The # threads to use. Default: 1	 # only for cnvkit

#### requires
- [`cnvkit`](http://cnvkit.readthedocs.io/en/stable/index.html)
- [`cnvnator`](https://github.com/abyzovlab/CNVnator)
- `wandy`: Inside cnv caller


###  pBam2FastqPE
#### description
- Convert sam/bam files to pair-end fastq files.

#### input
- `infile:file`: The sam/bam file. 
- - Sam files only available for biobambam, picard

#### output
- `fqfile1:file`: The 1st match of paired reads
- `fqfile2:file`: The 2nd match of paired reads

#### args
- `tool`                : The tool to use. Default: biobambam (bedtools, samtools, picard)
- `biobambam_bamtofastq`: The path of bamtofastq of biobambam. Default: bamtofastq
- `bedtools`            : The path of bedtools. Default: bedtools
- `samtools`            : The path of samtools. Default: samtools
- `picard`              : The path of picard. Default: picard
- `mem`                 : The memory to be used by picard. Default: 8G
- `gz`                  : Whether gzip the output files. Default: True
- `params`:             : Other params for `tool`. Default: ''
- `tmpdir`              : The tmpdir. Default: `__import__('tempfile').gettempdir()`

#### requires
- [picard](https://broadinstitute.github.io/picard/command-line-overview.html)
- [biobambam](https://github.com/gt1/biobambam2)
- [samtools](https://github.com/samtools/samtools)
- [bedtools](http://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html)


###  pBam2FastqSE
#### description
- Convert sam/bam files to single-end fastq files.

#### input
- `infile:file`: The sam/bam file. 
- - Sam files only available for biobambam, picard

#### output
- `fqfile:file`: The fastq file

#### args
- `tool`                : The tool to use. Default: biobambam (bedtools, samtools, picard)
- `biobambam_bamtofastq`: The path of bamtofastq of biobambam. Default: bamtofastq
- `bedtools`            : The path of bedtools. Default: bedtools
- `samtools`            : The path of samtools. Default: samtools
- `picard`              : The path of picard. Default: picard
- `mem`                 : The memory to be used by picard. Default: 8G
- `gz`                  : Whether gzip the output files. Default: True
- `params`:             : Other params for `tool`. Default: ''
- `tmpdir`              : The tmpdir. Default: `__import__('tempfile').gettempdir()`

#### requires
- [picard](https://broadinstitute.github.io/picard/command-line-overview.html)
- [biobambam](https://github.com/gt1/biobambam2)
- [samtools](https://github.com/samtools/samtools)
- [bedtools](http://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html)


###  pBam2Counts 

## GATK

###  pRealignerTargetCreator
#### description
- The local realignment process is designed to consume one or more BAM files and to locally realign reads such that the number of mismatching bases is minimized across all the reads. In general, a large percent of regions requiring local realignment are due to the presence of an insertion or deletion (indels) in the individual's genome with respect to the reference genome. Such alignment artifacts result in many bases mismatching the reference near the misalignment, which are easily mistaken as SNPs. Moreover, since read mapping algorithms operate on each read independently, it is impossible to place reads on the reference genome such that mismatches are minimized across all reads. Consequently, even when some reads are correctly mapped with indels, reads covering the indel near just the start or end of the read are often incorrectly mapped with respect the true indel, also requiring realignment. Local realignment serves to transform regions with misalignments due to indels into clean reads containing a consensus indel suitable for standard variant discovery approaches.
- Note that indel realignment is no longer necessary for variant discovery if you plan to use a variant caller that performs a haplotype assembly step, such as HaplotypeCaller or MuTect2. However it is still required when using legacy callers such as UnifiedGenotyper or the original MuTect. There are 2 steps to the realignment process:
- - Determining (small) suspicious intervals which are likely in need of realignment (RealignerTargetCreator)
- - Running the realigner over those intervals (see the IndelRealigner tool)
- For more details, see [the indel realignment method documentation](http://www.broadinstitute.org/gatk/guide/article?id=38).

#### input
- `bamfile:file`:  The aligned bam file
- `reffile`: The reference file

#### brings
- `bamfile`: `{{bamfile | bn}}.bai` The index file of input bam file
- `reffile#fai`: `{{reffile | bn}}.fai`
- `reffile#dict`: `{{reffile | bn}}.dict`

#### output
- `outfile:file`: A list of target intervals to pass to the IndelRealigner.

#### args
- `gatk`:     The gatk executable, default: "gatk"
- `picard`:   The picard executable, default: "picard"
- `params`:   Other parameters for RealignerTargetCreator, default: ""
- `samtools`: The samtools executable, default: "samtools"
- `tmpdir`:  The tmpdir to use. Default: /tmp
- `javamem`: The memory for java vm. Default: "-Xms1g -Xmx8g"

#### requires
- [GATK](https://software.broadinstitute.org/gatk)
- [samtools](http://www.htslib.org/) if `reffile` is not indexed or `bamfile` is not indexed.
- [picard](https://broadinstitute.github.io/picard/) if `reffile` is not dicted.


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
- `reffile:file`: The reference file

#### brings
- `bamfile`: `{{bamfile | bn}}.bai` The index file of input bam file
- `reffile#fai`: `{{reffile | bn}}.fai`
- `reffile#dict`: `{{reffile | bn}}.dict`

#### output
- `outfile:file`: A realigned version of input BAM file.

#### args
- `gatk`:     The gatk executable, default: "gatk"
- `picard`:   The picard executable, default: "picard"
- `params`:  Other parameters for IndelRealigner, default: ""
- `samtools`: The samtools executable, default: samtools
- `tmpdir`:  The tmpdir to use. Default: /tmp
- `javamem`: The memory for java vm. Default: "-Xms1g -Xmx8g"

#### requires
- [GATK](https://software.broadinstitute.org/gatk)
- [samtools](http://www.htslib.org/) if `reffile` is not indexed or `bamfile` is not indexed.
- [picard](https://broadinstitute.github.io/picard/) if `reffile` is not dicted.


###  pBaseRecalibrator  
#### description
- Variant calling algorithms rely heavily on the quality scores assigned to the individual base calls in each sequence read. These scores are per-base estimates of error emitted by the sequencing machines. Unfortunately the scores produced by the machines are subject to various sources of systematic technical error, leading to over- or under-estimated base quality scores in the data. Base quality score recalibration (BQSR) is a process in which we apply machine learning to model these errors empirically and adjust the quality scores accordingly. This allows us to get more accurate base qualities, which in turn improves the accuracy of our variant calls. The base recalibration process involves two key steps: first the program builds a model of covariation based on the data and a set of known variants (which you can bootstrap if there is none available for your organism), then it adjusts the base quality scores in the data based on the model. There is an optional but highly recommended step that involves building a second model and generating before/after plots to visualize the effects of the recalibration process. This is useful for quality control purposes. This tool performs the first step described above: it builds the model of covariation and produces the recalibration table. It operates only at sites that are not in dbSNP; we assume that all reference mismatches we see are therefore errors and indicative of poor base quality. This tool generates tables based on various user-specified covariates (such as read group, reported quality score, cycle, and context). Assuming we are working with a large amount of data, we can then calculate an empirical probability of error given the particular covariates seen at this site, where p(error) = num mismatches / num observations. The output file is a table (of the several covariate values, number of observations, number of mismatches, empirical quality score).

#### input
- `bamfile:file`: A BAM file containing data that needs to be recalibrated.
- `reffile:file`: The reference file

#### brings
- `bamfile`: `{{bamfile | bn}}.bai` The index file of input bam file
- `reffile#fai`: `{{reffile | bn}}.fai`
- `reffile#dict`: `{{reffile | bn}}.dict`

#### output
- `outfile:file`: A GATKReport file with many tables:
- - The list of arguments
- - The quantized qualities table
- - The recalibration table by read group
- - The recalibration table by quality score
- - The recalibration table for all the optional covariates

#### args
- `gatk`:     The gatk executable, default: "gatk"
- `params`:  Other parameters for BaseRecalibrator, default: ""
- `knownSites`: The known polymorphic sites to mask out, required
- `samtools`: The samtools executable, default: samtools
- `picard`:   The picard executable, default: "picard"
- `tmpdir`:  The tmpdir to use. Default: /tmp
- `javamem`: The memory for java vm. Default: "-Xms1g -Xmx8g"

#### requires
- [GATK](https://software.broadinstitute.org/gatk)
- [samtools](http://www.htslib.org/) if `reffile` is not indexed or `bamfile` is not indexed.
- [picard](https://broadinstitute.github.io/picard/) if `reffile` is not dicted.


###  pPrintReads   
#### description
- PrintReads is a generic utility tool for manipulating sequencing data in SAM/BAM format. It can dynamically merge the contents of multiple input BAM files, resulting in merged output sorted in coordinate order. It can also optionally filter reads based on various read properties such as read group tags using the `--read_filter/-rf` command line argument (see documentation on read filters for more information).
- Note that when PrintReads is used as part of the Base Quality Score Recalibration workflow, it takes the `--BQSR` engine argument, which is listed under Inherited Arguments > CommandLineGATK below.

#### input
- `bamfile:file`:    A BAM file.
- `recaltable:file`: The GATKReport file
- `reffile:file`:    The reference file

#### brings
- `bamfile`: `{{bamfile | bn}}.bai` The index file of input bam file
- `reffile#fai`: `{{reffile | bn}}.fai`
- `reffile#dict`: `{{reffile | bn}}.dict`

#### output
- `outfile:file`: A single processed bam file.

#### args
- `gatk`:     The gatk executable, default: "gatk"
- `params`:  Other parameters for PrintReads, default: ""
- `samtools`: The samtools executable, default: samtools
- `picard`:   The picard executable, default: "picard"
- `tmpdir`:  The tmpdir to use. Default: /tmp
- `javamem`: The memory for java vm. Default: "-Xms1g -Xmx8g"

#### requires
- [GATK](https://software.broadinstitute.org/gatk)
- [samtools](http://www.htslib.org/) if `reffile` is not indexed or `infile` is not indexed.
- [picard](https://broadinstitute.github.io/picard/) if `reffile` is not dicted.


###  pHaplotypeCaller 
#### description
- PrintReads is a generic utility tool for manipulating sequencing data in SAM/BAM format. It can dynamically merge the contents of multiple input BAM files, resulting in merged output sorted in coordinate order. It can also optionally filter reads based on various read properties such as read group tags using the `--read_filter/-rf` command line argument (see documentation on read filters for more information).
- Note that when PrintReads is used as part of the Base Quality Score Recalibration workflow, it takes the `--BQSR` engine argument, which is listed under Inherited Arguments > CommandLineGATK below.

#### input
- `bamfile:file`: A BAM file.
- `reffile:file`: The reference file

#### brings
- `bamfile`: `{{bamfile | bn}}.bai` The index file of input bam file
- `reffile#fai`: `{{reffile | bn}}.fai`
- `reffile#dict`: `{{reffile | bn}}.dict`

#### output
- `outfile:file`: Either a VCF or gVCF file with raw, unfiltered SNP and indel calls.

#### args
- `gatk`    : The gatk executable, default: "gatk"
- `params`  : Other parameters for HaplotypeCaller, default: ""
- `samtools`: The samtools executable, default: samtools
- `picard`:   The picard executable, default: "picard"
- `tmpdir`:  The tmpdir to use. Default: /tmp
- `javamem`: The memory for java vm. Default: "-Xms1g -Xmx8g"
- `nthread`: Corresponding to -nct option

#### requires
- [GATK](https://software.broadinstitute.org/gatk)
- [samtools](http://www.htslib.org/) if `reffile` is not indexed or `infile` is not indexed.
- [picard](https://broadinstitute.github.io/picard/) if `reffile` is not dicted.


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
- `reffile:file`: The reference file

#### brings
- `reffile#fai`: `{{reffile | bn}}.fai`
- `reffile#dict`: `{{reffile | bn}}.dict`

#### output
- `outfile:file`: A new VCF file containing the selected subset of variants.

#### args
- `gatk`:     The gatk executable, default: "gatk"
- `params`:  Other parameters for SelectVariants, default: ""
- `samtools`: The samtools executable, default: samtools
- `picard`:   The picard executable, default: "picard"
- `tmpdir`:  The tmpdir to use. Default: /tmp
- `javamem`: The memory for java vm. Default: "-Xms1g -Xmx8g"

#### requires
- [GATK](https://software.broadinstitute.org/gatk)
- [samtools](http://www.htslib.org/) if `reffile` is not indexed or `infile` is not indexed.
- [picard](https://broadinstitute.github.io/picard/) if `reffile` is not dicted.


###  pVariantFiltration
#### description
- This tool is designed for hard-filtering variant calls based on certain criteria. Records are hard-filtered by changing the value in the FILTER field to something other than PASS. Filtered records will be preserved in the output unless their removal is requested in the command line.
- The most common way of specifying filtering criteria is by using JEXL queries. See the article on JEXL expressions in the documentation Guide for detailed information and examples.

#### input
- `vcffile:file`: A variant call set from which to select a subset.
- `reffile:file`: The reference file

#### brings
- `reffile#fai`: `{{reffile | bn}}.fai`
- `reffile#dict`: `{{reffile | bn}}.dict`

#### output
- `outfile:file`: A filtered VCF.

#### args
- `gatk`:     The gatk executable, default: "gatk -T VariantFiltration"
- `params`:  Other parameters for VariantFiltration, default: ""
- `samtools`: The samtools executable, default: samtools
- `picard`:   The picard executable, default: "picard"
- `tmpdir`:  The tmpdir to use. Default: /tmp
- `javamem`: The memory for java vm. Default: "-Xms1g -Xmx8g"

#### requires
- [GATK](https://software.broadinstitute.org/gatk)
- [samtools](http://www.htslib.org/) if `reffile` is not indexed or `infile` is not indexed.
- [picard](https://broadinstitute.github.io/picard/) if `reffile` is not dicted.


###  pMuTect2
#### description
- MuTect2 is a somatic SNP and indel caller that combines the DREAM challenge-winning somatic genotyping engine of the original MuTect ([Cibulskis et al., 2013](http://www.nature.com/nbt/journal/v31/n3/full/nbt.2514.html)) with the assembly-based machinery of HaplotypeCaller. The basic operation of MuTect2 proceeds similarly to that of the HaplotypeCaller.
- NOTE: only Tumor/Normal variant calling implemented in bioprocs

#### input
- `tumor:file`:   the tumor bam file
- `normal:file`:  the normal bam file
- `reffile:file`: the reference file

#### brings
- `tumor`:  `{{tumor | bn}}.bai` the index file of tumor
- `normal`: `{{normal | bn}}.bai` the index file of normal
- `reffile#fai`: `{{reffile | bn}}.fai`
- `reffile#dict`: `{{reffile | bn}}.dict`

#### output
- `outfile:file`: The vcf file containing somatic mutations

#### args
- `gatk`:     The gatk executable, default: "gatk"
- `samtools`: The samtools executable, default: samtools
- `params`:   Other parameters for MuTect2, default: ""
- `picard`:   The picard executable, default: "picard"
- `tmpdir`:  The tmpdir to use. Default: /tmp
- `javamem`: The memory for java vm. Default: "-Xms1g -Xmx8g"

#### requires
- [GATK](https://software.broadinstitute.org/gatk)
- [samtools](http://www.htslib.org/) if index files of input files are not found
- [picard](https://broadinstitute.github.io/picard/) if `reffile` is not dicted.


###  pMuTect2Interval
#### description
- Use interval file model of MuTect2

#### input
- `tumor:file`:   the tumor bam file
- `normal:file`:  the normal bam file
- `reffile:file`: the reference file

#### brings
- `tumor`:  `{{tumor | bn}}.bai` the index file of tumor
- `normal`: `{{normal | bn}}.bai` the index file of normal
- `reffile#fai`: `{{reffile | bn}}.fai`
- `reffile#dict`: `{{reffile | bn}}.dict`

#### output
- `outfile:file`: The vcf file containing somatic mutations

#### args
- `gatk`:     The gatk executable, default: "gatk"
- `samtools`: The samtools executable, default: samtools
- `params`:   Other parameters for MuTect2, default: ""
- `picard`:   The picard executable, default: "picard"
- `tmpdir`:  The tmpdir to use. Default: /tmp
- `javamem`: The memory for java vm. Default: "-Xms1g -Xmx8g"

#### requires
- [GATK](https://software.broadinstitute.org/gatk)
- [samtools](http://www.htslib.org/) if index files of input files are not found
- [picard](https://broadinstitute.github.io/picard/) if `reffile` is not dicted.


## UTILS

## WXS

## COMMON

###  pSort
#### description
- Sort file using linux command `sort`

#### input
- `infile:file`: The input file

#### output
- `outfile:file`: The output file

#### args
- `skip`:   To skip first N lines. Default: 0
- `params`: The arguments used by `sort`


###  pFiles2Dir
#### description
- A helper process to convert a list of files into a directory, so that some processes can take it as input

#### input
- `infiles:files`: The input files

#### output
- `outdir:dir`:    The output directory


###  pFiles2List
#### description
- Put files to a list file

#### input
- `infiles:files`: The input files

#### args
- `delimit`: The delimit. Default: r"\n"

#### output
- `outfile:file`:  The output list file


###  pPat2Dir
#### description
- A helper process to convert a list of files by a pattern (wildcards) into a directory, so that some processes can take it as input

#### input
- `pattern:var`: The pattern

#### output
- `outdir:dir`:    The output directory


###  pMergeFiles
#### description
- Merge files in the input directory

#### input
- `indir:file`: The input directory

#### output
- `outfile:file`: The output file


###  pCbindList
#### description
- Column bind lists, fill miss rows with specific value

#### input
- `indir:file`: The directory containing the list files
- - header can be omited, but row names are required

#### output
- `outfile:file`: The output matrix

#### args
- `header`: Whether list has header. Default: False (will use file name as header)
- `na`:     The missing values. Default: 0
- - If it's a string, remember the quote (i.e.: '"missing"')


###  pFile2Proc
#### description
- Convert a file to a proc so it can be used as dependent

#### input
- `infile:file`: The input file

#### output
- `outfile:file`: The output file


## CNVKIT

###  pCNVkitAccess
#### description
- Calculate the sequence-accessible coordinates in chromosomes from the given reference genome, output as a BED file.

#### input
- `fafile:file`: The fasta file

#### output
- `outfile:file`: The output file

#### args
- `params`: Other parameters for `cnvkit.py access`
- `cnvkit`: The executable of cnvkit. Default: 'cnvkit.py'

#### requires
- [CNVkit](http://cnvkit.readthedocs.io/)


###  pCNVkitTarget
#### description
- Generate targets file for CNVkit using access file and annotate file (`cnvkit.py target`)

#### input
- `acfile:file`: The access file
- `anfile:file`: The annotate file

#### output
- `outfile:file`: The targets file

#### args
- `cnvkit`:  The executable of cnvkit. Default: 'cnvkit.py'
- `params`: Other parameters for `cnvkit.py target`

#### requires
- [CNVkit](http://cnvkit.readthedocs.io/)


###  pCNVkitCov
#### description
- Calculate coverage in the given regions from BAM read depths.

#### input
- `infile:file`: The bam file

#### output
- `outfile:file`: The output cnn file

#### args
- `tgfile`:  The target file
- `cnvkit`:  The executable of cnvkit. Default: 'cnvkit.py'
- `nthread`: The number of threads to use. Default: 1
- `params`:  Other parameters for `cnvkit.py coverage`

#### requires
- [CNVkit](http://cnvkit.readthedocs.io/)


###  pCNVkitRef
#### description
- Compile a copy-number reference from the given files or directory (containing normal samples). If given a reference genome (-f option), also calculate the GC content and repeat-masked proportion of each region.

#### input
- `indir:file`:  The input directory containing the cnn files

#### output
- `outfile:file`: The output reference cnn file

#### args
- `cnvkit`:  The executable of cnvkit. Default: 'cnvkit.py'
- `params`:  Other parameters for `cnvkit.py reference`, default: " --no-edge "

#### requires
- [CNVkit](http://cnvkit.readthedocs.io/)


###  pCNVkitFix
#### description
- Combine the uncorrected target and antitarget coverage tables (.cnn) and correct for biases in regional coverage and GC content, according to the given reference. Output a table of copy number ratios (.cnr)

#### input
- `infile:file`:  The cnn file to be fixed
- `rcfile:file`:  The reference cnn file

#### output
- `outfile:file`: The cnr file

#### args
- `cnvkit`:  The executable of cnvkit. Default: 'cnvkit.py'
- `params`:  Other parameters for `cnvkit.py fix`, default: " --no-edge "

#### requires
- [CNVkit](http://cnvkit.readthedocs.io/)


###  pCNVkitSeg
#### description
- Infer discrete copy number segments from the given coverage table

#### input
- `infile:file`:  The cnr file 

#### output
- `outfile:file`: The cns file

#### args
- `cnvkit`:  The executable of cnvkit. Default: 'cnvkit.py'
- `nthread`: The number of threads to use. Default: 1
- `params`:  Other parameters for `cnvkit.py segment`, default: ""

#### requires
- [CNVkit](http://cnvkit.readthedocs.io/)


###  pCNVkitCall
#### description
- Given segmented log2 ratio estimates (.cns), derive each segment's absolute integer copy number 

#### input
- `infile:file`:  The cns file 

#### output
- `outfile:file`: The callcns file

#### args
- `cnvkit`:  The executable of cnvkit. Default: 'cnvkit.py'
- `params`:  Other parameters for `cnvkit.py segment`, default: ""

#### requires
- [CNVkit](http://cnvkit.readthedocs.io/)


###  pCNVkitPlot
#### description
- Plot CNVkit results

#### input
- `cnrdir:file`:  The directory containing copy number ratio files
- `cnsdir:file`:  The directory containing copy number segment files

#### output
- `outdir:dir`:   The output directory

#### args
- `cnvkit`:   The executable of cnvkit. Default: 'cnvkit.py'
- `region`:       The region for zoom-in plots. Default: '' (don't plot zoom-in view)
- `gene`:         The genes to be highlighted. Default: ''
- `scatter`:      Whether to generate the scatter plot. Default: True
- `diagram`:      Whether to generate the diagram plot. Default: True
- `heatmap`:      Whether to generate the heatmap plot. Default: True

#### requires
- [CNVkit](http://cnvkit.readthedocs.io/)


###  pCNVkitRpt
#### description
- Report CNVkit results

#### input
- `cnrfile:file`:  The file containing copy number ratio
- `cnsfile:file`:  The file containing copy number segment

#### output
- `outdir:dir`:   The output directory

#### args
- `cnvkit`:   The executable of cnvkit. Default: 'cnvkit.py'
- `breaks`:       Whether to report breakpoints. Default: True
- `gainloss`:     Whether to report gainloss. Default: True
- `metrics`:      Whether to report metrics. Default: True
- `segmetrics`:   Whether to report segmetrics. Default: True

#### requires
- [CNVkit](http://cnvkit.readthedocs.io/)


###  pCNVkit2Vcf
#### description
- Output vcf file for cnvkit results

#### input
- `cnsfile:file`: The cns file

#### output
- `outfile:file`: The vcf file

#### args
- `cnvkit`:   The executable of cnvkit. Default: 'cnvkit.py'
- `params`:   Other params for `cnvkit.py export`

#### requires
- [CNVkit](http://cnvkit.readthedocs.io/)


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


###  pGetPromoterBed
#### description
- Get the promoter region in bed format

#### input
- `gene`: the gene

#### output
- `outfile:file`: the bed file containing the promoter region

#### args
- `up`: the upstream to the tss, default: 2000
- `down`: the downstream to the tss, default: 2000
- `genome`: the genome, default: hg19

#### require
- [python-mygene](http://mygene.info/)


###  pGetPromotersBed
#### description
- Get the promoter regions in bed format of a gene list give in genefile

#### input
- `genefile:file`: the gene list file

#### output
- `outfile:file`: the bed file containing the promoter region

#### args
- `up`: the upstream to the tss, default: 2000
- `down`: the downstream to the tss, default: 2000
- `genome`: the genome, default: hg19

#### require
- [python-mygene](http://mygene.info/)


## TCGA

###  pDownload
#### description
- Download TCGA use `gdc-client` and a manifest file

#### input
- `manifile:file`: the manifest file

#### output
- `outdir:file`:   the directory containing downloaded file

#### args
- `params`:        other params for `gdc-client`, default: "--no-related-files --no-file-md5sum -n 20"
- `bin-gdc`:       the executable file of `gdc-client`, default: "gdc-client"


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


## WXSSTAT

###  pVcf2List
#### description
- Convert vcf to stat files for pCallRate

#### input
- `vcffile:file`: The vcf file

#### output
- `outfile:file`: The stat file

#### args
- `chroms`: SNPs on chromosomes to consider, default: "" (all chroms)
- - use "chr1-22, chrX, chrY" for chr1 to chr22, chrX and chrY

#### requires
- [`pyvcf`](https://github.com/jamescasbon/PyVCF)


###  pCallRate
#### description
- Calculate sample/snp call rate from a matrix of snp-sample
- - rows are snps, columns are samples

#### input
- `infile:file`:    The snp-sample matrix file

#### output
- `outsample:file`: The report of call rate for each sample
- `figsample:file`: The bar chat of sample call rates
- `outsnp:file`:    The report of call rate for each snp
- `figsnp:file`:    The bar chat of snp call rates


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


###  pSnpEff2Stat
#### description
- Convert csvstat file from snpEff to R-readable matrix for plotting

#### input
- `indir:file`: The directory containing the csv stat files from `snpEff ann`

#### output
- `outdir:dir`: The output directory

#### args
- `chroms`:     The chromsome filter. Default: "" (all chroms)
- - Note: snpEff csvstat file has no "chr" prefix


###  pPlotSnpEff
#### description
- Plot snpEff annotation statistics

#### input
- `indir:file`: The snpEff result directory containing matrix files generated by pSnpEff2Stat

#### output
- `outdir:dir`: The output directory

#### requires
- [`pwwang/corrplot`](https://github.com/pwwang/corrplot)
- - use `library(devtools); install.github("pwwang/corrplot")`
- [`ggplot2`](http://ggplot2.org/)


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

