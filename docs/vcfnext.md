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

## pCallRate

### description
Calculate sample/snp call rate from single sample vcfs

### input
#### `indir:file`::     The dir containing the vcfs  

## pCepip

### input
#### `infile:file`:: The input file (vcf or avinput)  

### output
#### `outfile:file`:: The cepip result file  

### args
#### `cepip`::    The path of cepip  
#### `cell` ::    The related cell line  
#### `params`::   Other params for cepip  

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

## pMafMerge

### description
Merge maf files.

### input
#### `indir:dir`:: The directory containing the maf files  

### output
#### `outfile:file`:: The merged maf file  

## pMaftools

### args
#### `ngenes`::   
{% endraw %}
