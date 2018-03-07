# tcga
<!-- toc -->
{% raw %}

## pDownload

### description
Download TCGA use `gdc-client` and a manifest file

### input
#### `manifile:file`:: the manifest file  

### output
#### `outdir:file`::   the directory containing downloaded file  

## pSample2SubmitterID

### description
convert TCGA sample names with submitter id with metadata and sample containing folder

### input
#### `dir:file`::    the directory containing the samples  
#### `mdfile:file`:: the metadata file  

## pConvertExpFiles2Matrix

### description
convert TCGA expression files to expression matrix, and convert sample name to submitter id

### input
#### `dir:file`::    the directory containing the samples  
#### `mdfile:file`:: the metadata file  

### output
#### `outfile:file`::the output matrix  

## pConvertMutFiles2Matrix

### description
convert TCGA mutation files (vcf.gz) to mut matrix, and convert sample name to submitter id

### input
#### `dir:file`::    the directory containing the samples  
#### `mdfile:file`:: the metadata file  
{% endraw %}
