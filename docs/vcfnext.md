{% raw %}

## pVcfStatsPlot

### description
	Convert csvstat file from snpEff to R-readable matrix and plot them.

### input
#### `indir:file`:
 The directory containing the csv stat files from `snpEff ann`  

### output
#### `outdir:dir`:
 The output directory  

## pCallRate

### description
	Calculate sample/snp call rate from single sample vcfs

### input
#### `indir:file`:
     The dir containing the vcfs  

## pCepip

### input
#### `infile:file`:
 The input file (vcf or avinput)  

### output
#### `outfile:file`:
 The cepip result file  

### args
#### `cepip`:
    The path of cepip  
#### `cell` :
    The related cell line  
#### `params`:
   Other params for cepip  
{% endraw %}
