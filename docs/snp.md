
## pSnp2Bed

### description
	Find coordinates for SNPs in BED format.

### input
#### `snpfile:file`:
 the snp file, each snp per line  

### output
#### `outfile:file`:
 the result file, columns are:  
		- chrom, start(0-based), end, name, score, strand, ref, allele

### args
#### `genome`:
 default: hg19  
#### `snpver`:
 default: snp147  

## pSnp2Avinput

### description
	Convert SNP list to avinput to ANNOVAR.

### input
#### `snpfile:file`:
 the snp file, each snp per line  

### output
#### `outfile:file`:
 the result avinput file  

### args
#### `genome`:
 default: hg19  
#### `snpver`:
 default: snp147  
