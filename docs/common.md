{% raw %}

## pSort

### description
	Sort file using linux command `sort`

### input
#### `infile:file`:
 The input file  

### output
#### `outfile:file`:
 The output file  

## pFiles2Dir

### description
	A helper process to convert a list of files into a directory, so that some processes can take it as input

### input
#### `infiles:files`:
 The input files  

## pFile2Proc

### description
	Convert a file to a proc so it can be used as dependent

### input
#### `infile:file`:
 The input file  

## pStr2File

### description
	Save string to a file.

### input
#### `in:var`:
 The input string.  

## pAddHeader

### description
	Add the header of 1st file to 2nd file.

### input
#### `infile1:file`:
 The first file containing the header.  
#### `infile2:file`:
 The second file with the body.  

### output
#### `outfile:file`:
 The output file with the header from 1st input file, body from 2nd file.  
{% endraw %}
