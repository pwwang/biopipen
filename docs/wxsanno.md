<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [pSnpEff](#psnpeff)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->


## pSnpEff

### description
	This is the default command. It is used for annotating variant filed (e.g. VCF files).

### input
#### `infile:file`:
  The input file   

### output
#### `outdir:file`:
 The directory containing output anntated file, snpEff_genes.txt and snpEff_summary.html  

### args
#### `snpEff`:
       The snpEff executable, default: "snpEff"  
#### `params`:
    Other parameters for `snpEff`, default: "-Xms1g -Xmx4g -v"  
#### `genome`:
    The genome used for annotation, default: "hg19"  
#### `informat`:
  The format of input file [vcf or bed], default: "vcf"  
#### `outformat`:
 The format of output file [vcf, gatk, bed, bedAnn], default: "vcf"  
#### `csvStats`:
  Whether to generate csv stats file, default: True.  
#### `htmlStats`:
 Whether to generate the html summary file, default: False.  
#### `javamem`:
   The memory to use. Default: '-Xms1g -Xmx8g'  
