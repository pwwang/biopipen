# snp
<!-- toc -->
{% raw %}

## pSnp2Bedx

### description
Find coordinates for SNPs in BEDX format.

### input
#### `snpfile:file`:: the snp file, each snp per line  

### output
#### `outfile:file`:: the result file, columns are:  
	- chrom, start(0-based), end, name, score, strand, ref, allele

### args
#### `genome`:: default: hg19  
#### `snpver`:: default: snp147  
#### `notfound`:: What to do if the snp is not found. Default: skip  
#### `inmeta`:: The metadata for input file to determine which column is rsID  
#### `xcols`:: The extra columns to extract and output to extra columns in output file.  
#### `indem`:: The input delimit. Default: '\\t'  
#### `incom`:: The input comment. Default: '#'  
#### `skip`:: The lines to skip for input file. Default: 0  

### requires
[`python-cruzdb`](https://github.com/brentp/cruzdb)

## pSnp2Avinput

### description
Convert SNP list to avinput to ANNOVAR.

### input
#### `snpfile:file`:: the snp file, each snp per line  

### output
#### `outfile:file`:: the result avinput file  

### args
#### `genome`:: default: hg19  
#### `snpver`:: default: snp147  

### requires
[`python-cruzdb`](https://github.com/brentp/cruzdb)
{% endraw %}
