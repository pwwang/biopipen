# cnvkit
<!-- toc -->
{% raw %}

## pCNVkitAccess

### description
Calculate the sequence-accessible coordinates in chromosomes from the given reference genome, output as a BED file.

### input
#### `fafile:file`:: The fasta file  

### output
#### `outfile:file`:: The output file  

### args
#### `params`:: Other parameters for `cnvkit.py access`  
#### `cnvkit`:: The executable of cnvkit. Default: 'cnvkit.py'  

## pCNVkitTarget

### description
Generate targets file for CNVkit using access file and annotate file (`cnvkit.py target`)

### input
#### `acfile:file`:: The access file  
#### `anfile:file`:: The annotate file  

### output
#### `outfile:file`:: The targets file  

### args
#### `cnvkit`::  The executable of cnvkit. Default: 'cnvkit.py'  
#### `params`:: Other parameters for `cnvkit.py target`  

## pCNVkitCov

### description
Calculate coverage in the given regions from BAM read depths.

### input
#### `infile:file`:: The bam file  

### output
#### `outfile:file`:: The output cnn file  

### args
#### `tgfile`::  The target file  
#### `cnvkit`::  The executable of cnvkit. Default: 'cnvkit.py'  
#### `nthread`:: The number of threads to use. Default: 1  
#### `params`::  Other parameters for `cnvkit.py coverage`  

## pCNVkitRef

### description
Compile a copy-number reference from the given files or directory (containing normal samples). If given a reference genome (-f option), also calculate the GC content and repeat-masked proportion of each region.

### input
#### `indir:file`::  The input directory containing the cnn files  

### output
#### `outfile:file`:: The output reference cnn file  

### args
#### `cnvkit`::  The executable of cnvkit. Default: 'cnvkit.py'  
#### `params`::  Other parameters for `cnvkit.py reference`, default: " --no-edge "  

## pCNVkitFix

### description
Combine the uncorrected target and antitarget coverage tables (.cnn) and correct for biases in regional coverage and GC content, according to the given reference. Output a table of copy number ratios (.cnr)

### input
#### `infile:file`::  The cnn file to be fixed  
#### `rcfile:file`::  The reference cnn file  

### output
#### `outfile:file`:: The cnr file  

### args
#### `cnvkit`::  The executable of cnvkit. Default: 'cnvkit.py'  
#### `params`::  Other parameters for `cnvkit.py fix`, default: " --no-edge "  

## pCNVkitSeg

### description
Infer discrete copy number segments from the given coverage table

### input
#### `infile:file`::  The cnr file   

### output
#### `outfile:file`:: The cns file  

### args
#### `cnvkit`::  The executable of cnvkit. Default: 'cnvkit.py'  
#### `nthread`:: The number of threads to use. Default: 1  
#### `params`::  Other parameters for `cnvkit.py segment`, default: ""  

## pCNVkitCall

### description
Given segmented log2 ratio estimates (.cns), derive each segment's absolute integer copy number 

### input
#### `infile:file`::  The cns file   

### output
#### `outfile:file`:: The callcns file  

### args
#### `cnvkit`::  The executable of cnvkit. Default: 'cnvkit.py'  
#### `params`::  Other parameters for `cnvkit.py segment`, default: ""  

## pCNVkitPlot

### description
Plot CNVkit results

### input
#### `cnrdir:file`::  The directory containing copy number ratio files  
#### `cnsdir:file`::  The directory containing copy number segment files  

### output
#### `outdir:dir`::   The output directory  

### args
#### `cnvkit`::   The executable of cnvkit. Default: 'cnvkit.py'  
#### `region`::       The region for zoom-in plots. Default: '' (don't plot zoom-in view)  
#### `gene`::         The genes to be highlighted. Default: ''  
#### `scatter`::      Whether to generate the scatter plot. Default: True  
#### `diagram`::      Whether to generate the diagram plot. Default: True  
#### `heatmap`::      Whether to generate the heatmap plot. Default: True  

## pCNVkitRpt

### description
Report CNVkit results

### input
#### `cnrfile:file`::  The file containing copy number ratio  
#### `cnsfile:file`::  The file containing copy number segment  

### output
#### `outdir:dir`::   The output directory  

### args
#### `cnvkit`::   The executable of cnvkit. Default: 'cnvkit.py'  
#### `breaks`::       Whether to report breakpoints. Default: True  
#### `gainloss`::     Whether to report gainloss. Default: True  
#### `metrics`::      Whether to report metrics. Default: True  
#### `segmetrics`::   Whether to report segmetrics. Default: True  

## pCNVkit2Vcf

### description
Output vcf file for cnvkit results

### input
#### `cnsfile:file`:: The cns file  

### output
#### `outfile:file`:: The vcf file  

### args
#### `cnvkit`::   The executable of cnvkit. Default: 'cnvkit.py'  
#### `params`::   Other params for `cnvkit.py export`  
{% endraw %}
