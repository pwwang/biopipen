# fastx
<!-- toc -->
{% raw %}

## pFastqSim

### description
Simulate reads

### input
#### `seed`:
The seed to generate simulation file  
	- None: use current timestamp.

### output
#### `fq1:file`:
The first pair read file  
#### `fq2:file`:
The second pair read file  

### args
#### `tool`:
The tool used for simulation. Default: wgsim (dwgsim)  
#### `len1`:
The length of first pair read. Default: 100  
#### `len2`:
The length of second pair read. Default: 100  
#### `num`:
The number of read PAIRs. Default: 1000000  
#### `gz`:
Whether generate gzipped read file. Default: True  
#### `wgsim`:
The path of wgsim. Default: wgsim  
#### `dwgsim`:
The path of wgsim. Default: dwgsim  
#### `ref`:
The reference genome. Required  
#### `params`:
Other params for `tool`. Default: ""  

## pFastQC

### description
QC report for fastq file

### input
#### `fq:file`:
The fastq file (also fine with gzipped)  

### output
#### `outdir:dir`:
The output direcotry  

### args
#### `tool`:
The tool used for simulation. Default: fastqc   
#### `fastqc`:
The path of fastqc. Default: fastqc  
#### `nthread`:
Number of threads to use. Default: 1  
#### `params`:
Other params for `tool`. Default: ""  

## pFastMC

### description
Multi-QC based on pFastQC

### input
#### `qcdir:file`:
The direcotry containing QC files  

### output
#### `outdir:dir`:
The output direcotry  

### args
#### `tool`:
The tool used for simulation. Default: multiqc   
#### `multiqc`:
The path of fastqc. Default: multiqc  
#### `params`:
Other params for `tool`. Default: ""  

## pFastqTrim

### description
Trim pair-end FASTQ reads

### input
#### `fq1:file`:
The input fastq file  
#### `fq2:file`:
The input fastq file  

### output
#### `outfq1:file`:
The trimmed fastq file  
#### `outfq2:file`:
The trimmed fastq file  

### args
#### `tool`        :
The tools used for trimming. Default: trimmomatic (cutadapt|skewer)  
#### `cutadapt`    :
The path of seqtk. Default: cutadapt  
#### `skewer`      :
The path of fastx toolkit trimmer. Default: skewer  
#### `trimmomatic` :
The path of trimmomatic. Default: trimmomatic  
#### `params`      :
Other params for `tool`. Default: ""  
#### `nthread`     :
Number of threads to be used. Default: 1  
- Not for cutadapt
#### `gz`          :
Whether gzip output files. Default: True  
#### `mem`         :
The memory to be used. Default: 4G  
- Only for trimmomatic
#### `minlen`      :
Discard trimmed reads that are shorter than `minlen`. Default: 18  
- For trimmomatic, the number will be `minlen`*2 for MINLEN, as it filters before trimming
#### `minq`        :
Minimal mean qulity for 4-base window or leading/tailing reads. Default: 3  
#### `cut5`        :
Remove the 5'end reads if they are below qulity. Default: 3  
#### `cut3`        :
Remove the 3'end reads if they are below qulity. Default: 3  
- Not for skewer
#### `adapter1`    :
The adapter for sequence. Default: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC  
#### `adapter2`    :
The adapter for pair-end sequence. Default: AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA  

## pFastqSETrim

### description
Trim single-end FASTQ reads

### input
#### `fq:file`:
The input fastq file  

### output
#### `outfq:file`:
The trimmed fastq file  

### args
#### `tool`        :
The tools used for trimming. Default: trimmomatic (cutadapt|skewer)  
#### `cutadapt`    :
The path of seqtk. Default: cutadapt  
#### `skewer`      :
The path of fastx toolkit trimmer. Default: skewer  
#### `trimmomatic` :
The path of trimmomatic. Default: trimmomatic  
#### `params`      :
Other params for `tool`. Default: ""  
#### `nthread`     :
Number of threads to be used. Default: 1  
- Not for cutadapt
#### `gz`          :
Whether gzip output files. Default: True  
#### `mem`         :
The memory to be used. Default: 4G  
- Only for trimmomatic
#### `minlen`      :
Discard trimmed reads that are shorter than `minlen`. Default: 18  
- For trimmomatic, the number will be `minlen`*2 for MINLEN, as it filters before trimming
#### `minq`        :
Minimal mean qulity for 4-base window or leading/tailing reads. Default: 3  
#### `cut5`        :
Remove the 5'end reads if they are below qulity. Default: 3  
#### `cut3`        :
Remove the 3'end reads if they are below qulity. Default: 3  
- Not for skewer
#### `adapter`     :
The adapter for sequence. Default: AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC  

## pFastqSE2Sam

### description
Cleaned paired fastq (.fq, .fq.gz, .fastq, .fastq.gz file to mapped sam/bam file

## pFastq2Sam

### description
Cleaned paired fastq (.fq, .fq.gz, .fastq, .fastq.gz file to mapped sam/bam file
{% endraw %}
