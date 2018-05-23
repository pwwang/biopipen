# common
<!-- toc -->
{% raw %}

## pSort

### description
Sort file using linux command `sort`

### input
#### `infile:file`:: The input file  

### output
#### `outfile:file`:: The output file  

### args
#### `skip`::   To skip first N lines. Default: 0  
#### `case`::   Case-sensitivity. Default: True  
	- If True, will set $LANG as C
	- Otherwise, $LANG will be set as en_US.UTF-8
#### `mem`    :: The buffer size. Default: 4G  
#### `tmpdir` :: The tmpdir.  
#### `unique` :: Just keep the unique lines. Default: False  
#### `delimit`:: The delimit to separate the fields. Default: '\t'  
#### `params` :: The arguments used by `sort`  

## pFiles2Dir

### description
A helper process to convert a list of files into a directory, so that some processes can take it as input

### input
#### `infiles:files`:: The input files  

### output
#### `outdir:dir`::    The output directory  

## pFile2Proc

### description
Convert a file to a proc so it can be used as dependent

### input
#### `infile:file`:: The input file  

### output
#### `outfile:file`:: The output file  

## pStr2File

### description
Save string to a file.

### input
#### `in:var`:: The input string.  

### output
#### `outfile:file`:: The output file.  

## pHead

### description
Get the top N lines from a file

### input
#### `infile:file`:: The input file  

### output
#### `outfile:file`:: The output file  

### args
#### `n`:: Top n lines. You may use '-n' to skip last n lines.  

## pTail

### description
Get the bottom N lines from a file

### input
#### `infile:file`:: The input file  

### output
#### `outfile:file`:: The output file  

### args
#### `n`:: Bottom n lines. You may use '+n' to skip first n lines.  

## pPrepend

### description
Prepend a string to a file

### input
#### `in:var`:: The input string.  
#### `infile:file`:: The input file.  

### output
#### `outfile:file`:: The output file.  

## pAddHeader

### description
Add the header of 1st file to 2nd file.

### input
#### `infile1:file`:: The first file containing the header.  
#### `infile2:file`:: The second file with the body.  

### output
#### `outfile:file`:: The output file with the header from 1st input file, body from 2nd file.  

### args
#### `n`:: The number of header lines.  

## pMergeFiles

### description
Merge files in the input directory

### input
#### `indir:file`:: The input directory  

### output
#### `outfile:file`:: The output file  

### args
#### `inopts`:: The options for input file.  
	- defaults: skip: 0, comment: #, delimit '\\t'
#### `outopts`:: The options for output file. Defaults:  
	- head: False (not output head line)
	- headPrefix: `#` (The prefix for head line)
	- headDelimit: `\\t` (The delimit for head line)
	- headTransform: `None` (The callback for head line)
	- delimit: `\\t` (The delimit for output line)

## pSplitRows

### description
Split a file by rows, specially usefull to split a job into multithreads/multiprocesses.

### input
#### `infile:file`:: The input file  

### output
#### `outdir:dir`:: The output directory including the split files  

### args
#### `skip`:: The skip first n lines. Default: `0`  
#### `cnames`:: The column names. If True, the column names will be added to each split file. Default: `True`  
#### `n`:: Number of files to split. Default: `8`  
{% endraw %}
