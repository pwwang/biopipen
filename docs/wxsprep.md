{% raw %}

## pTrimmomaticPE

### description
	Trimming Illumina NGS paired-end data

### input
#### `fqfile1:file`:
 The 1st fastq file (could be in .gz format)  
#### `fqfile2:file`:
 The 2nd fastq file  

### output
#### `outfile1:file`:
 The 1st output file  
#### `outfile2:file`:
 The 2nd output file  

### args
#### `trimmomatic`:
    The trimmomatic executable, default: "trimmomatic"  
#### `phred`:
  "phred33" (default) or "phred64"  
#### `params`:
 Other params for trimmomatric, default: "ILLUMINACLIP:{adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"  
		- have to replace `{adapter}` with the path of the adapter file
#### `nthread`:
 1  

## pTrimmomaticSE

### description
	Trimming Illumina NGS single-end data

### input
#### `fqfile:file`:
 The fastq file (could be in .gz format)  

### output
#### `outfile:file`:
 The output file  

### args
#### `trimmomatic`:
    The trimmomatic executable, default: "trimmomatic"  
#### `phred`:
  "phred33" (default) or "phred64"  
#### `params`:
 Other params for trimmomatric, default: "ILLUMINACLIP:{adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36"  
		- have to replace `{adapter}` with the path of the adapter file
#### `nthread`:
 1  

## pAlignPEByBWA

### description
	Align paired-end reads to reference genome using bwa mem

### input
#### `infile1:file`:
 read file 1 (fastq, or fastq gzipped)  
#### `infile2:file`:
 read file 2 (fastq, or fastq gzipped)  
#### `reffile:file`:
 The reference file  

### output
#### `outfile:file`:
 The output sam file  

### args
#### `bwa`:
    The bwa executable, default: bwa  
#### `params`:
 Other params for bwa mem, default: "-M"  
#### `nthread`:
 1  

## pAlignSEByBWA

### description
	Align paired-end reads to reference genome using bwa mem

### input
#### `infile:file`:
  read file (fastq, or fastq gzipped)  
#### `reffile:file`:
 The reference file  

### brings
#### `reffile#bwt`:
 "{{reffile | bn}}.bwt",   
#### `reffile#sa`:
  "{{reffile | bn}}.sa",  
#### `reffile#ann`:
 "{{reffile | bn}}.ann",  
#### `reffile#amb`:
 "{{reffile | bn}}.amb",  
#### `reffile#pac`:
 "{{reffile | bn}}.pac"  

### output
#### `outfile:file`:
 The output sam file  

### args
#### `bwa`:
    The bwa executable, default: bwa  
#### `params`:
 Other params for bwa mem, default: "-M"  
#### `nthread`:
 1  
#### `reffile`:
 The reference file, required  

## pAlignPEByNGM

### description
	Align paired-end reads to reference genome using NextGenMap

### input
#### `infile1:file`:
 read file 1 (fastq, or fastq gzipped)  
#### `infile2:file`:
 read file 2 (fastq, or fastq gzipped)  
#### `reffile:file`:
 The reference file  

### output
#### `outfile:file`:
 The output sam/bam file  

### args
#### `ngm`:
    The NextGenMap executable, default: ngm  
#### `nthread`:
 1  
#### `outtype`:
 sam or bam, default: sam (only sam for now, due to bug of ngm 0.5.3 (fixed in 0.5.4))  
#### `params`:
 Other params for ngm, default: "--rg-id ngm --rg-sm sample"  

## pAlignSEByNGM

### description
	Align single-end reads to reference genome using NextGenMap

### input
#### `infile1:file`:
 read file 1 (fastq, or fastq gzipped)  
#### `infile2:file`:
 read file 2 (fastq, or fastq gzipped)  
#### `reffile:file`:
 The reference file  

### output
#### `outfile:file`:
 The output sam/bam file  

### args
#### `ngm`:
    The NextGenMap executable, default: ngm  
#### `nthread`:
 1  
#### `outtype`:
 sam or bam, default: sam (only sam for now, due to bug of ngm 0.5.3 (fixed in 0.5.4))  
#### `params`:
 Other params for ngm, default: "--rg-id ngm --rg-sm sample"  

## pMergeBams

### description
	Merge bam files

### input
#### `bamdir:dir`:
   the dir containing bam files   

### output
#### `outfile:file`:
 the merged bam file  

### args
#### `samtools`:
 the executable path of samtools, default: "samtools"  
#### `nthread`:
      Number of BAM/CRAM compression threads  
#### `params`:
       Other parameters for `samtools merge`, default: ""  
{% endraw %}
