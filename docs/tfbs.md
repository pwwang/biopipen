<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [pMotifScan](#pmotifscan)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->


## pMotifScan

### input
#### `tffile:file`:
 The infile containing TF name and motif name.  
		- If only one column is give, will be used as both TF and motif name
		- If there are 2 columns, 1st column will be motif name, 2nd column will be TF name
#### `sfile:file`:
 The sequence file  

### output
#### `outdir:file`:
 The output dir  

### args
#### `tools`   :
 The tool used to scan the motif. Default: 'meme'  
#### `meme`    :
 The path of MEME's fimo. Default: 'fimo'  
	`motifs   : The motif database in MEME format.
#### `pval`    :
 The pvalue cutoff. Default: 1e-4  
#### `cleanmname`:
 Whether to clean motif name. Default: True  
#### `ucsclink`:
 The ucsc link template. Default: `https://genome.ucsc.edu/cgi-bin/hgTracks?db=hg19&position={}`  
#### `nthread` :
 Number of threads used to scan, only available when you have multiple mids. Default: 1  
#### `params`  :
 Other parameters for `fimo`  
