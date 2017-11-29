<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [pSnp2Bed](#psnp2bed)
- [pSnp2Avinput](#psnp2avinput)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->


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
