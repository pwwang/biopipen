# matrix
<!-- toc -->
{% raw %}

## pMatrix

### description
Operate a matrix and save the new matrix to file.

### input
#### `infile:file`:
The input file containing the matrix  

### output
#### `outfile:file`:
The output matrix  

## pCbind

### description
Cbind the rest of files to the first file.

### input
#### `infiles:files`:
The input files  

### output
#### `outfile:file`:
The output matrix  

## pRbind

### description
Rbind the rest of files to the first file.

### input
#### `infiles:files`:
The input files  

### output
#### `outfile:file`:
The output matrix  

## pCsplit

### description
Split a matrix by columns and save them into files.

### input
#### `infile:file`:
The input file  

### output
#### `outdir:dir`:
The directory containing the output column files  

## pRsplit

### description
Split a matrix by rows and save them into files.

### input
#### `infile:file`:
The input file  

### output
#### `outdir:dir`:
The directory containing the output row files  

## pTxtFilter

### description
Filter a tab-delimit file (txt/tsv file)

### input
#### `infile:file`:
The input file  

### output
#### `outfile:file`:
The output file  

## pTxtTransform

### description
Transform a tab-delimit file (txt/tsv file)

### input
#### `infile:file`:
The input file  

### output
#### `outfile:file`:
The output file  

## pSimRead

### description
Read files simultaneously.
NOTE: only one file allows multiple lines with same value to compare, and that file should be the first one. For example: 
```
File1:
1	1
1	2
1	3
File2:
1	1
2	2
3	3
```
If you compare the first column, File1 has to put at the begining for input.

### input
#### `infiles:files`:
The input files  

### output
#### `outfile:file`:
The output file  

### args
#### `skip`:
argument skip for each file  
#### `delimit`:
argument delimit for each file  
#### `usehead`:
The header from which input file will be used for output file.  
	- Default: None (Don't write header)
#### `gzip`:
argument gzip for each file  
#### `match`:
The match function.   
#### `do`:
The do function. Global vaiable `fout` is available to write results to output file.  
{% endraw %}
