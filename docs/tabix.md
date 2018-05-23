# tabix
<!-- toc -->
{% raw %}

## pTabix

### description
Use tabix to extract information.

### input
#### `infile`:: a local or remote file  
#### `region`:: a region or a file containing regions  

### output
#### `outfile:file`:: The information extracted from the input file  

### args
#### `tabix`:: The path to `tabix`  
#### `params`:: Other params for `tabix`  

## pTabixIndex

### description
Generate tabix index file.

### input
#### `infile:file`:: the input file  
	- Could be bgzipped.

### output
#### `outfile:file`:: The bgzipped file  
#### `outidx:file`::  The tabix index file  

### args
#### `tabix`:: The path to `tabix`  
#### `params`:: Other params for `tabix`  
#### `python`:: Will be used to generate command line arguments.  
{% endraw %}
