<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [pCNVnator](#pcnvnator)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->


## pCNVnator

### description
	Use `CNVnator` to call CNVs from bam file

### input
#### `infile:file`:
  The bam file   

### output
#### `outfile:file`:
 The vcf file  

### args
#### `cnvnator`:
      The CNVnator executable, default: "cnvnator"  
#### `cnv2vcf`:
  The converter executable to convert CNVnator results to vcf, default: "cnvnator2VCF.pl"  
#### `binsize`:
  The bin_size, default: 100  
#### `genome`:
   The genome: default: hg19  
#### `chrom`:
    Chromosome names, default: "" (all chromosomes)  
#### `chrdir`:
   The dir contains reference sequence of chromosomes, default: "" (don't specify)  
	
