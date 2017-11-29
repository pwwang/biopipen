<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [pVcfStatsPlot](#pvcfstatsplot)
- [pCallRate](#pcallrate)
- [pCepip](#pcepip)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->


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
