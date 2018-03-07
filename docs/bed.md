# bed
<!-- toc -->
{% raw %}

## pBedSort

### description
Sort bed files

### input
#### `infile:file`:: The input file  

### output
#### `outfile:file`:: The output file  

### args
#### `tool`::         The tool used to sort the file. Default: sort (bedtools, bedops)  
#### `bedtools`::     The path to bedtools. Default: bedtools  
#### `bedops_sort`::  The path to bedops' sort-bed. Default: sort-bed  
#### `mem`::          The memory to use. Default: 8G  
#### `tmpdir`::       The tmpdir to use. Default: `$TMPDIR`  
#### `unique`::       Remove the dupliated records? Default: True  
#### `params`::       Other params for `tool`. Default: {}  

## pBedCluster

### description
Assign cluster id to each record

### input
#### `infile:file`:: The input bed file  

### output
#### `outfile:file`:: The output file  

### args
#### `tool`::         The tool used to sort the file. Default: bedtools  
#### `bedtools`::     The path to bedtools. Default: bedtools  
#### `params`::       Other params for `tool`. Default: ''  
{% endraw %}
