<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [pDownload](#pdownload)
- [pSample2SubmitterID](#psample2submitterid)
- [pConvertExpFiles2Matrix](#pconvertexpfiles2matrix)
- [pConvertMutFiles2Matrix](#pconvertmutfiles2matrix)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->


## pDownload

### description
	Download TCGA use `gdc-client` and a manifest file

### input
#### `manifile:file`:
 the manifest file  

### output
#### `outdir:file`:
   the directory containing downloaded file  

## pSample2SubmitterID

### description
	convert TCGA sample names with submitter id with metadata and sample containing folder

### input
#### `dir:file`:
    the directory containing the samples  
#### `mdfile:file`:
 the metadata file  

## pConvertExpFiles2Matrix

### description
	convert TCGA expression files to expression matrix, and convert sample name to submitter id

### input
#### `dir:file`:
    the directory containing the samples  
#### `mdfile:file`:
 the metadata file  

### output
#### `outfile:file`:
the output matrix  

## pConvertMutFiles2Matrix

### description
	convert TCGA mutation files (vcf.gz) to mut matrix, and convert sample name to submitter id

### input
#### `dir:file`:
    the directory containing the samples  
#### `mdfile:file`:
 the metadata file  
