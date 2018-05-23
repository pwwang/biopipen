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

### args
#### `params`::        other params for `gdc-client`, default: "--no-related-files --no-file-md5sum -n 20"  
#### `bin-gdc`::       the executable file of `gdc-client`, default: "gdc-client"  

## pSample2SubmitterID

### description
convert TCGA sample names with submitter id with metadata and sample containing folder

### input
#### `dir:file`::    the directory containing the samples  
#### `mdfile:file`:: the metadata file  

### output
#### `outdir:file`:: the directory containing submitter-id named files  

## pConvertExpFiles2Matrix

### description
convert TCGA expression files to expression matrix, and convert sample name to submitter id

### input
#### `dir:file`::    the directory containing the samples  
#### `mdfile:file`:: the metadata file  

### output
#### `outfile:file`::the output matrix  

### requires
[python-mygene](https://pypi.python.org/pypi/mygene/3.0.0)

## pConvertMutFiles2Matrix

### description
convert TCGA mutation files (vcf.gz) to mut matrix, and convert sample name to submitter id

### input
#### `dir:file`::    the directory containing the samples  
#### `mdfile:file`:: the metadata file  

### output
#### `outfile:file`::the output matrix  
{% endraw %}
