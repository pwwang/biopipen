# vcfnext
<!-- toc -->
{% raw %}

## pVcfStatsPlot

### description
Convert csvstat file from snpEff to R-readable matrix and plot them.

### input
#### `indir:file`:: The directory containing the csv stat files from `snpEff ann`  

### output
#### `outdir:dir`:: The output directory  

### args
#### `chroms`::     The chromsome filter. Default: "" (all chroms)  
- Note: snpEff csvstat file has no "chr" prefix

## pCallRate

### description
Calculate sample/snp call rate from single sample vcfs

### input
#### `indir:file`::     The dir containing the vcfs  

### output
#### `outsample:file`:: The report of call rate for each sample  
#### `figsample:file`:: The bar chat of sample call rates  
#### `outsnp:file`::    The report of call rate for each snp  
#### `figsnp:file`::    The bar chat of snp call rates  

## pCepip

### input
#### `infile:file`:: The input file (vcf or avinput)  

### output
#### `outfile:file`:: The cepip result file  

### args
#### `cepip`::    The path of cepip  
#### `cell` ::    The related cell line  
#### `params`::   Other params for cepip  

### requires
[`cepip`](http://jjwanglab.org/cepip/)

## pMutSig

### description
MutSig stands for "Mutation Significance".  MutSig analyzes lists of mutations discovered in DNA sequencing, to identify genes that were mutated more often than expected by chance given background mutation processes.
For more information, see Lawrence, M. et al. Mutational heterogeneity in cancer and the search for new cancer-associated genes. Nature 499, 214-218 (2013).

See [dcumentation](http://archive.broadinstitute.org/cancer/cga/mutsig_run)

### input
#### `infile:file`:: mutation table  

### output
#### `outdir:dir`:: The output directory  

### args
#### `mutsig` :: The path to `run_MutSigCV.sh`, default: 'mutsig'  
#### `mcr`    :: The Matlab MCR path  
#### `cvrg`   :: coverage table  
#### `cvrt`   :: covariates table  
#### `mutdict`:: mutation_type_dictionary_file  
#### `chrdir` :: chr_files_hg18 or chr_files_hg19  

### requires
[MutSig](http://archive.broadinstitute.org/cancer/cga/mutsig_download)

## pMafMerge

### description
Merge maf files.

### input
#### `indir:dir`:: The directory containing the maf files  

### output
#### `outfile:file`:: The merged maf file  

### args
#### `excols`:: How to deal with extra columns other than 34 standard columns from TCGA.  
	- merge(default): Merge the columns, if one not exists, fill with an empty string.
	- discard: Just discard the extra columns, with only 34 columns left. So you can also put just one maf file in the indir with some columns missed to fill it with standard columns.

## pMaftools

### args
#### `ngenes`::   

### requires
[``]
{% endraw %}
