# fastx
<!-- toc -->
{% raw %}

## pFastqSim

### description
Simulate reads

### input
#### `seed`:: The seed to generate simulation file  
	- None: use current timestamp.

### output
#### `fq1:file`:: The first pair read file  
#### `fq2:file`:: The second pair read file  

### args
#### `tool`::  The tool used for simulation. Default: wgsim (dwgsim)  
#### `len1`::  The length of first pair read. Default: 100  
#### `len2`::  The length of second pair read. Default: 100  
#### `num`::   The number of read PAIRs. Default: 1000000  
#### `gz`::    Whether generate gzipped read file. Default: True  
#### `wgsim`:: The path of wgsim. Default: wgsim  
#### `dwgsim`::The path of wgsim. Default: dwgsim  
#### `ref`::   The reference genome. Required  
#### `params`::Other params for `tool`. Default: ""  

### requires
[`wgsim`](https://github.com/lh3/wgsim)

## pFastQC

### description
QC report for fastq file

### input
#### `fq:file`::    The fastq file (also fine with gzipped)  

### output
#### `outdir:dir`:: The output direcotry  

### args
#### `tool`::    The tool used for simulation. Default: fastqc  
#### `fastqc`::  The path of fastqc. Default: fastqc  
#### `nthread`:: Number of threads to use. Default: 1  
#### `params`::Other params for `tool`. Default: ""  

### requires
[`fastqc`](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

## pFastMC

### description
Multi-QC based on pFastQC

### input
#### `qcdir:file`::  The direcotry containing QC files  

### output
#### `outdir:dir`:: The output direcotry  

### args
#### `tool`::    The tool used for simulation. Default: multiqc  
#### `multiqc`:: The path of fastqc. Default: multiqc  
#### `params`::  Other params for `tool`. Default: ""  

### requires
[`multiqc`](http://multiqc.info/)

## pFastqTrim

### description
Trim pair-end FASTQ reads

### input
#### `fq1:file`::  The input fastq file  
#### `fq2:file`::  The input fastq file  

### output
#### `outfq1:file`:: The trimmed fastq file  
#### `outfq2:file`:: The trimmed fastq file  

### args
#### `tool`        :: The tools used for trimming. Default: trimmomatic (cutadapt|skewer)  
#### `cutadapt`    :: The path of seqtk. Default: cutadapt  
#### `skewer`      :: The path of fastx toolkit trimmer. Default: skewer  
#### `trimmomatic` :: The path of trimmomatic. Default: trimmomatic  
#### `params`      :: Other params for `tool`. Default: ""  
#### `nthread`     :: Number of threads to be used. Default: 1  
- Not for cutadapt
#### `gz`          :: Whether gzip output files. Default: True  
#### `mem`         :: The memory to be used. Default: 4G  
- Only for trimmomatic
#### `minlen`      :: Discard trimmed reads that are shorter than `minlen`. Default: 18  
- For trimmomatic, the number will be `minlen`*2 for MINLEN, as it filters before trimming
#### `minq`        :: Minimal mean qulity for 4-base window or leading/tailing reads. Default: 3  
#### `cut5`        :: Remove the 5'end reads if they are below qulity. Default: 3  
#### `cut3`        :: Remove the 3'end reads if they are below qulity. Default: 3  
- Not for skewer
#### `adapter1`    :: The adapter for sequence. Default: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC  
#### `adapter2`    :: The adapter for pair-end sequence. Default: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA  

### requires
[`cutadapt`](http://cutadapt.readthedocs.io/en/stable/guide.html)
[`skewer`](https://github.com/relipmoc/skewer)
[`trimmomatic`](https://github.com/timflutre/trimmomatic)

## pFastqSETrim

### description
Trim single-end FASTQ reads

### input
#### `fq:file`::  The input fastq file  

### output
#### `outfq:file`:: The trimmed fastq file  

### args
#### `tool`        :: The tools used for trimming. Default: trimmomatic (cutadapt|skewer)  
#### `cutadapt`    :: The path of seqtk. Default: cutadapt  
#### `skewer`      :: The path of fastx toolkit trimmer. Default: skewer  
#### `trimmomatic` :: The path of trimmomatic. Default: trimmomatic  
#### `params`      :: Other params for `tool`. Default: ""  
#### `nthread`     :: Number of threads to be used. Default: 1  
- Not for cutadapt
#### `gz`          :: Whether gzip output files. Default: True  
#### `mem`         :: The memory to be used. Default: 4G  
- Only for trimmomatic
#### `minlen`      :: Discard trimmed reads that are shorter than `minlen`. Default: 18  
- For trimmomatic, the number will be `minlen`*2 for MINLEN, as it filters before trimming
#### `minq`        :: Minimal mean qulity for 4-base window or leading/tailing reads. Default: 3  
#### `cut5`        :: Remove the 5'end reads if they are below qulity. Default: 3  
#### `cut3`        :: Remove the 3'end reads if they are below qulity. Default: 3  
- Not for skewer
#### `adapter`     :: The adapter for sequence. Default: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC  

### requires
[`cutadapt`](http://cutadapt.readthedocs.io/en/stable/guide.html)
[`skewer`](https://github.com/relipmoc/skewer)
[`trimmomatic`](https://github.com/timflutre/trimmomatic)

## pFastqSE2Sam

### description
Cleaned paired fastq (.fq, .fq.gz, .fastq, .fastq.gz file to mapped sam/bam file

### args
#### `tool`::   The tool used for alignment. Default: bwa (bowtie2|ngm)  
#### `bwa`::    Path of bwa, default: bwa  
#### `ngm`::    Path of ngm, default: ngm  
#### `bowtie2`::Path of bowtie2, default: bowtie2  
#### `rg`::     The read group. Default: {'id': '', 'pl': 'Illumina', 'pu': 'unit1', 'lb': 'lib1', 'sm': ''}  
- `id` will be parsed from filename with "_LX_" in it if not given
- `sm` will be parsed from filename
#### `ref`::    Path of reference file  
#### `params`:: Other params for tool, default: ''  

## pFastq2Sam

### description
Cleaned paired fastq (.fq, .fq.gz, .fastq, .fastq.gz file to mapped sam/bam file

### args
#### `tool`   :: The tool used for alignment. Default: bwa (bowtie2, ngm, star)  
#### `bwa`    :: Path of bwa, default: bwa  
#### `ngm`    :: Path of ngm, default: ngm  
#### `star`   :: Path of ngm, default: STAR  
#### `bowtie2`:: Path of bowtie2, default: bowtie2  
#### `rg`::     The read group. Default: {'id': '', 'pl': 'Illumina', 'pu': 'unit1', 'lb': 'lib1', 'sm': ''}  
- `id` will be parsed from filename with "_LX_" in it if not given
- `sm` will be parsed from filename
#### `ref`    :: Path of reference file  
#### `refgene`:: The GTF file for STAR to build index. It's not neccessary if index is already been built. Default: ''  
#### `params` :: Other params for tool, default: ''  
{% endraw %}
