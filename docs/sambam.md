# sambam
<!-- toc -->
{% raw %}

## pSam2Bam

### description
Deal with mapped sam/bam files, including sort, markdup, and/or index

### input
#### `infile:file`:: The input file  

### output
#### `outfile:file`:: The output bam file  
#### `idxfile:file`:: The index of the output bam file  
- If args.index == False, it'll a link to outfile and should be never used

### args
#### `tool`             :: The tool used to do the sort. Default: sambamba (picard|sambamba|biobambam|samtools)  
#### `sambamba`         :: The path of the sambamba. Default: sambamba  
#### `picard`           :: The path of the picard. Default: picard  
#### `biobambam_bamsort`:: The path of the biobambam's bamsort. Default: bamsort  
#### `samtools`         :: The path of the samtools. Default: samtools  
#### `sort`             :: Do sorting? Default: True  
- If input is sam, tool is biobambam, this should be True
#### `index`            :: Do indexing? Default: True  
#### `markdup`          :: Do duplicates marking? Default: False  
- `rmdup` for samtools will be called
#### `rmdup`            :: Do duplicates removing? Default: False  
#### `tmpdir`           :: The tmp dir used to store tmp files. Default: <system default tmpdir>  
#### `sortby`           :: Sort by coordinate or queryname. Default: coordinate  
#### `nthread`          :: Default: 1  
#### `informat`         :: The format of input file. Default: <detect from extension> (sam|bam)  
#### `params`           :: Other parameters for `tool`. Defaut: ""  
#### `mem`              :: The max memory to use. Default: "16G"  
- Unit could be G/g/M/m
- Will be converted to -Xmx4G, and -Xms will be 1/8 of it

### requires
[sambamba](https://lomereiter.github.io/sambamba/docs/sambamba-view.html) if `args.tool` == samtools or reference used but not indexed.
[picard](https://broadinstitute.github.io/picard/command-line-overview.html)
[biobambam](https://github.com/gt1/biobambam2)
[samtools](https://github.com/samtools/samtools)

## pBamMarkdup

### description
Mark/remove duplicates for bam files

### input
#### `infile:file`:: The input file  

### output
#### `outfile:file`:: The output bam file  

### args
#### `tool`             :: The tool used to do the sort. Default: sambamba (picard|sambamba|biobambam|samtools|bamutil)  
#### `sambamba`         :: The path of sambamba. Default: sambamba  
#### `picard`           :: The path of picard. Default: picard  
#### `biobambam_bamsort`:: The path of biobambam's bamsort. Default: bamsort  
#### `samtools`         :: The path of samtools. Default: samtools  
#### `bamutil`          :: The path of bamutil. Default: bam  
#### `rmdup`            :: Do duplicates removing? Default: False  
- Samtools will anyway remove the duplicates
#### `tmpdir`           :: The tmp dir used to store tmp files. Default: <system default tmpdir>  
#### `nthread`          :: Default: 1  
- Not available for samtools and picard
#### `params`           :: Other parameters for `tool`. Defaut: ""  
#### `mem`              :: The max memory to use. Default: "16G"  
- Unit could be G/g/M/m
- Will be converted to -Xmx4G, and -Xms will be 1/8 of it

### requires
[sambamba](https://lomereiter.github.io/sambamba/docs/sambamba-view.html)
[picard](https://broadinstitute.github.io/picard/command-line-overview.html)
[biobambam](https://github.com/gt1/biobambam2)
[samtools](https://github.com/samtools/samtools)
[bamutil](http://genome.sph.umich.edu/wiki/BamUtil#Programs)

## pBamRecal

### description
Recalibrate a bam file

### input
#### `infile:file`:: The bam file  

### output
#### `outfile:file`:: The output bam file  

### args
#### `tool`                         :: The tool used to recalibrate the bam file. Default: `gatk` (gatk|bamutil)  
#### `gatk`                         :: The path of gatk, including java path. Default: `gatk`  
#### `samtools`                     :: The path of samtools. Default: `samtools`  
#### `bamutil`                      :: The path of bamutil. Default: `bam`  
#### `picard`                       :: The path of picard. Default: `picard`  
#### `paramsRealignerTargetCreator` :: Other parameters for `gatk RealignerTargetCreator`. Defaut: ""  
#### `paramsIndelRealigner`         :: Other parameters for `gatk IndelRealigner`. Defaut: ""  
#### `paramsBaseRecalibrator`       :: Other parameters for `gatk BaseRecalibrator`. Defaut: ""  
#### `paramsPrintReads`             :: Other parameters for `gatk PrintReads`. Defaut: ""  
#### `params`                       :: Other parameters for `bam recab`. Default: ""  
#### `mem`                          :: The max memory to use. Default: "32G"  
#### `knownSites`                   :: The known polymorphic sites to mask out. Default: "" (Required for GATK)  
#### `ref`                          :: The reference file. Required.  
- Will be converted to -Xmx4G, and -Xms will be 1/8 of it

### requires
[gatk](https://software.broadinstitute.org/gatk)
[samtools](https://github.com/samtools/samtools) if `args.ref` is not indexed, or bamutil is used for bam index file generation.
[picard](https://broadinstitute.github.io/picard/command-line-overview.html) if `args.ref is not dicted.`

## pBamReadGroup

### description
Add or replace read groups of a bam file

### input
#### `infile:file`:: The bam file  

### output
#### `outfile:file`:: The output bam file  

### args
#### `tool`                         :: The tool used. Default: `picard` (picard|bamutil)  
#### `picard`                       :: The path of picard. Default: `picard`  
#### `bamutil`                      :: The path of bamutil. Default: `bam`  
#### `rg`                           :: The read group. Default: {'id': '', 'pl': 'Illumina', 'pu': 'unit1', 'lb': 'lib1', 'sm': ''}  
- `id` will be parsed from filename with "_LX_" in it if not given
- `sm` will be parsed from filename
#### `params`                       :: Other parameters for `tool`. Defaut: ""  
#### `mem`                          :: The max memory to use. Default: "4G"  
- Will be converted to -Xmx4G, and -Xms will be 1/8 of it
#### `tmpdir`                       :: The temporary directory. Default: <system tmpdir>  

### requires
[gatk](https://lomereiter.github.io/sambamba/docs/sambamba-view.html)
[samtools](https://github.com/samtools/samtools) if `args.ref` is not indexed.
[picard](https://broadinstitute.github.io/picard/command-line-overview.html) if `args.ref is not dicted.`

## pBamReorder

### description
Reorder a sam/bam file by a given reference file using `picard ReorderSam`

### input
#### `infile:file`:: The sam/bam file  

### output
#### `outfile:file`:: The output bam file  

### args
#### `picard`                       :: The path of picard. Default: `picard`  
#### `ref`                          :: The reference file. Required  
#### `params`                       :: Other parameters for `picard ReorderSam`. Defaut: ""  
#### `mem`                          :: The max memory to use. Default: "4G"  
- Will be converted to -Xmx4G, and -Xms will be 1/8 of it
#### `tmpdir`                       :: The temporary directory. Default: <system tmpdir>  

### requires
[picard](https://broadinstitute.github.io/picard/command-line-overview.html)

## pBamMerge

### description
Merges multiple SAM and/or BAM files (must be sorted by coordinate) into a single file.

### input
#### `infiles:file`:: Input sam/bam files to be merged  

### output
#### `outfile:file`:: The merged bam file  

### args
#### `tool`     :: The tool used to merge. Default: bamutil (picard|samtools|sambamba)  
#### `picard`   :: The path of picard. Default: `picard`  
#### `bamutil`  :: The path of bamutil. Default: `bam`  
#### `samtools` :: The path of samtools. Default: `samtools`  
#### `sambamba` :: The path of sambamba. Default: `sambamba`  
#### `params`   :: Other parameters for `tool`. Defaut: ""  
#### `mem`      :: The max memory to use. Default: "4G"  
- Will be converted to -Xmx4G, and -Xms will be 1/8 of it, just for picard
#### `tmpdir`   :: The temporary directory. Default: <system tmpdir>  
#### `nthread`  :: # threads to use. Default: 1  
- For picard, if nthread>1, USE_THREADING=true, otherwise USE_THREADING=false

### requires
[picard](https://broadinstitute.github.io/picard/command-line-overview.html)

## pBam2Gmut

### description
Call germline (snps and indels) from a call-ready bam file.

### input
#### `infile:file`:: The input bam file  

### output
#### `outfile:file`:: The vcf file containing the mutations  

### args
#### `tool`::         The tool used to call mutations. Default: gatk (vardict, snvsniffer, platypus, strelka)  
#### `gatk`::         The path of gatk. Default: gatk  
#### `vardict`::      The path of vardict. Default: vardict  
#### `snvsniffer`::   The path of snvsniffer. Default: SNVSniffer  
#### `samtools`::     The path of samtools. Default: samtools (used to generate reference index)  
#### `platypus`::     The path of platypus. Default: platypus  
#### `strelka`::      The path of strelka. Default: configureStrelkaGermlineWorkflow.py  
#### `configParams`:: The params for `strelka` configuration. Default: ""  
#### `picard`::       The path of picard. Default: picard  
#### `mem`::          The memory to be used. Default: 32G  
- will be converted to -Xms4G -Xmx32G for java programs
#### `ref`::          The reference file. Required.  
#### `gz`::           Gzip output file? Default: False  
#### `tmpdir`::       The temporary directory. Default: <system tmpdir>  
#### `params`::       Other params for `tool`. Default: ""  

### requires
[gatk](https://lomereiter.github.io/sambamba/docs/sambamba-view.html)
[samtools](https://github.com/samtools/samtools) if `args.ref` is not indexed.
[picard](https://broadinstitute.github.io/picard/command-line-overview.html) if `args.ref is not dicted.`
[vardict](https://github.com/AstraZeneca-NGS/VarDict)
[snvsniffer](http://snvsniffer.sourceforge.net/homepage.htm#latest)
[platypus](http://www.well.ox.ac.uk/platypus)
[strelka@2.7.1+](https://github.com/Illumina/strelka)

## pBamPair2Smut

### description
Call somatic mutations from tumor-normal bam pair.

### input
#### `tumor:file`:: The tumor bam file  
#### `normal:file`:: The normal bam file  

### output
#### `outfile:file`:: The vcf file  

### args
#### `tool`:: The tool used to call mutations. Default: gatk (somaticsniper, strelka, snvsniffer, virmid, varidct)  
#### `gatk`:: The path to gatk. Default: gatk  
#### `somaticsniper`:: The path to gatk. Default: bam-somaticsniper  
#### `strelka`:: The path to gatk. Default: configureStrelkaSomaticWorkflow.py  
#### `snvsniffer`:: The path to gatk. Default: SNVSniffer  
#### `virmid`:: The path to gatk. Default: virmid  
#### `vardict`:: The path to gatk. Default: vardict  
#### `samtools`:: The path to gatk. Default: samtools  
#### `picard`:: The path to gatk. Default: picard  
#### `configParams`:: The configuration parameters for `configureStrelkaSomaticWorkflow.py`. Default: `{}`  
#### `params`:: The parameters for main programs. Default: `{}`  
#### `meme`:: The memory. Default: 24G  
#### `ref`:: The reference genom. Default: `params.ref.value`  
#### `gz`:: Whether gzip the output vcf file. Default: False  
#### `nthread`:: The number of threads to use. Default: 1  
#### `tmpdir`:: The temporary directory. Default: `params.tmpdir.value`  

### requires
[gatk](https://lomereiter.github.io/sambamba/docs/sambamba-view.html)
[samtools](https://github.com/samtools/samtools) if `args.ref` is not indexed.
[picard](https://broadinstitute.github.io/picard/command-line-overview.html) if `args.ref is not dicted.`
[vardict](https://github.com/AstraZeneca-NGS/VarDict)
[snvsniffer](http://snvsniffer.sourceforge.net/homepage.htm#latest)
[platypus](http://www.well.ox.ac.uk/platypus)
[strelka@2.7.1+](https://github.com/Illumina/strelka)

## pBam2Cnv

### description
Detect copy number variation from bam files.

### input
#### `input:file`:: The bam file  

### output
#### `outfile:file`:: The output vcf file  
#### `outdir`:: The output directory containing other result files  

### args
#### `gz`                    :: Whether to gzip the output vcf file. Default: False  
#### `tool`                  :: The tool used to call cnv. Default: 'cnvkit'  
#### `cnvnator`              :: The path of cnvnator. Default: 'cnvnator'  
#### `cnvnator2vcf`          :: The path of cnvnator2VCF. Default: 'cnvnator2VCF.pl'  
#### `cnvkit`                :: The path of cnvkit. Default: 'cnvkit.py'  
#### `wandy`                 :: Tha path of Wandy. Default: 'Wandy'. A `tool.info` file should be with the executable file.  
#### `ref`                   :: The reference file. Required by cnvkit to generate access file. Default: ''  
#### `cnvkitAccessParams`    :: The params for cnvkit access command. Default: '-s 5000'  
#### `cnvkitTargetParams`    :: The params for cnvkit target command. Default: '--split --short-names'  
#### `cnvkitCoverageParams`  :: The params for cnvkit coverage command. Default: ''  
#### `cnvkitReferenceParams` :: The params for cnvkit reference command. Default: '--no-edge'  
#### `cnvkitFixParams`       :: The params for cnvkit fix command. Default: '--no-edge'  
#### `cnvkitSegmentParams`   :: The params for cnvkit segment command. Default: ''  
#### `cnvkitCallParams`      :: The params for cnvkit call command. Default: ''  
#### `cnvkitPlotParams`      :: The params for cnvkit plot command. Default: ''  
#### `cnvkitBreaksParams`    :: The params for cnvkit breaks command. Default: ''  
#### `cnvkitGainlossParams`  :: The params for cnvkit gainloss command. Default: ''  
#### `cnvkitMetricsParams`   :: The params for cnvkit metrics command. Default: ''  
#### `cnvkitSegmetricsParams`:: The params for cnvkit segmetrics command. Default: '--iqr'  
#### `cnvkitExportParams`    :: The params for cnvkit export command. Default: ''  
#### `cnvkitScatterParams`   :: The params for cnvkit scatter command. Default: [''] # multiple scatter plots  
#### `cnvkitHeatmapParams`   :: The params for cnvkit heatmap command. Default: [''] # multiple heatmap plots  
#### `cnvkitDiagramParams`   :: The params for cnvkit diagram command. Default: ''  
#### `cnvkitReport`          :: Generate cnvkit reports? Default: True  
#### `cnvkitPlot`            :: Generate cnvkit plots? Default: True  
#### `cnvnatorBinsize`       :: Bin size for cnvnator. Default: 100  
#### `cnvnatorGenome`        :: Genome for cnvnator. Default: 'hg19'. (NCBI36, hg18, GRCh37, hg19)  
#### `params`                :: The params for `tool`. Default: '-t 1' # wandy 1:hg19 solid cell/blood, 2:hg19 cell free/plamsa, 3:hg38 solid cell/blood, 4:hg38 cell free/plamsa  
#### `mem`                   :: The memory used. Default: '20G' # only for wandy  
#### `nthread`               :: The # threads to use. Default: 1	 # only for cnvkit  

### requires
[`cnvkit`](http://cnvkit.readthedocs.io/en/stable/index.html)
[`cnvnator`](https://github.com/abyzovlab/CNVnator)
#### `wandy`:: Inside cnv caller  

## pBamStats

### description
Get read depth from bam files.

### input
#### `infile:file`:: The input bam file  

### output
#### `outfile:file`:: The output statistic file  
#### `outdir:dir`::   The directory containing result files and figures.  

### args
#### `tool`:: The tool used to do the job. Default: bamstats  
#### `bamstats`:: The path to bamstats. Default: bamstats  
#### `params`:: Other params to main program. Default: `{}`  
#### `mem`:: The memory to be used. Default: 16G  
#### `plot`:: Whether plot the result. Default: True  

## pBam2Fastq

### description
Convert sam/bam files to pair-end fastq files.

### input
#### `infile:file`:: The sam/bam file.  
	- Sam files only available for biobambam, picard

### output
#### `fqfile1:file`:: The 1st match of paired reads  
#### `fqfile2:file`:: The 2nd match of paired reads  

### args
#### `tool`     :: The tool to use. Default: biobambam (bedtools, samtools, picard)  
#### `biobambam`:: The path of bamtofastq of biobambam. Default: bamtofastq  
#### `bedtools` :: The path of bedtools. Default: bedtools  
#### `samtools` :: The path of samtools. Default: samtools  
#### `picard`   :: The path of picard. Default: picard  
#### `mem`      :: The memory to be used by picard. Default: 8G  
#### `gz`       :: Whether gzip the output files. Default: True  
#### `params`::  : Other params for `tool`. Default: ''  
#### `tmpdir`   :: The tmpdir. Default: `__import__('tempfile').gettempdir()`  

### requires
[picard](https://broadinstitute.github.io/picard/command-line-overview.html)
[biobambam](https://github.com/gt1/biobambam2)
[samtools](https://github.com/samtools/samtools)
[bedtools](http://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html)

## pBam2FastqSE

### description
Convert sam/bam files to single-end fastq files.

### input
#### `infile:file`:: The sam/bam file.  
	- Sam files only available for biobambam, picard

### output
#### `fqfile:file`:: The fastq file  

### args
#### `tool`     :: The tool to use. Default: biobambam (bedtools, samtools, picard)  
#### `biobambam`:: The path of bamtofastq of biobambam. Default: bamtofastq  
#### `bedtools` :: The path of bedtools. Default: bedtools  
#### `samtools` :: The path of samtools. Default: samtools  
#### `picard`   :: The path of picard. Default: picard  
#### `mem`      :: The memory to be used by picard. Default: 8G  
#### `gz`       :: Whether gzip the output files. Default: True  
#### `params`::  : Other params for `tool`. Default: ''  
#### `tmpdir`   :: The tmpdir. Default: `__import__('tempfile').gettempdir()`  

### requires
[picard](https://broadinstitute.github.io/picard/command-line-overview.html)
[biobambam](https://github.com/gt1/biobambam2)
[samtools](https://github.com/samtools/samtools)
[bedtools](http://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html)

## pBam2Counts

### description
Extract read counts from RNA-seq bam files.

### input
#### `infile:file`:: The input bam files  

### outfile
#### `outfile:file`:: The count file  

### args
#### `tool`:: The tool used to extract counts. Default: ht-seq  
#### `htseq`:: The path of htseq-count.  
#### `params`:: Other params for main program.  
#### `refgene`:: The reference gene in GTF format.  

### requires
[`htseq`](https://htseq.readthedocs.io/)
{% endraw %}
