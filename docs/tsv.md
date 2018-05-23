# tsv
<!-- toc -->
{% raw %}

## pMatrixR

### description
Operate a matrix and save the new matrix to file.

### input
#### `infile:file`:: The input file containing the matrix  

### output
#### `outfile:file`:: The output matrix  

### args
#### `cnames`:: Whether the input file has cnames. Default: True  
#### `rnames  `:: Whether the input file has rnames  . Default: 1  
#### `code`:: The R code to operating the matrix. (the matrix is read in variable `mat`)  

## pCbind

### description
Cbind the rest of files to the first file.

### input
#### `infiles:files`:: The input files  

### output
#### `outfile:file`:: The output matrix  

### args
#### `cnames`:: Whether the input file has cnames. Default: True  
	- or [True, True, False] corresponding to the file order
#### `rnames  `:: Whether the input file has rnames  . Default: 1  
#### `miss`:: Replacement for missing values. Default: `NA`  

## pRbind

### description
Rbind the rest of files to the first file.

### input
#### `infiles:files`:: The input files  

### output
#### `outfile:file`:: The output matrix  

### args
#### `cnames`:: Whether the input file has cnames. Default: True  
	- or [True, True, False] corresponding to the file order
#### `rnames  `:: Whether the input file has rnames  . Default: 1  
#### `miss`:: Replacement for missing values. Default: `NA`  

## pCsplit

### description
Split a matrix by columns and save them into files.

### input
#### `infile:file`:: The input file  

### output
#### `outdir:dir`:: The directory containing the output column files  

### args
#### `cnames`:: Whether the input file has cnames. Default: True  
	- or [True, True, False] corresponding to the file order
#### `rnames  `:: Whether the input file has rnames  . Default: 1  

## pRsplit

### description
Split a matrix by rows and save them into files.

### input
#### `infile:file`:: The input file  

### output
#### `outdir:dir`:: The directory containing the output row files  

### args
#### `cnames`:: Whether the input file has cnames. Default: True  
	- or [True, True, False] corresponding to the file order
#### `rnames  `:: Whether the input file has rnames  . Default: 1  

## pTsv

### description
Read, Transform, filter a TSV file.

### input
#### `infile:file`:: The input file  

### output
#### `outfile:file`:: The output file  

### args
#### `inopts`:: The input options for infile:  
	- `delimit`: The delimit. Default: `\\t`
	- `comment`: The comment sign. Default: `#`
	- `skip`: First N lines to skip. Default: `0`
	- `ftype`: The file type. Metadata can be assigned direct (list/OrderedDict). If not specified, metadata will be generated automatically.
#### `outopts`:: The output options for outfile:  
	- `delimit`: The delimit for records. Default: `\\t`
	- `head`: Output header or not. Default: `False`
	- `headDelimit`: The delimit for header. Default: `\\t`
	- `headPrefix`: The prefix for header. Default: ``
	- `headTransform`: The transformer for header. Default: `None`
	- `ftype`: The file type. Metadata can be assigned direct (list/OrderedDict, '+' as an element or key is allowed to indicate extra meta from the reader). If not specified, metadata will be borrowed from the reader. 
#### `ops`:: A ops function to transform the row. Argument is an instance of `readRecord`  
#### `opshelper`:: A helper function for `args.ops`  

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
#### `infiles:files`:: The input files  

### output
#### `outfile:file`:: The output file  

### args
#### `skip`:: argument skip for each file  
#### `delimit`:: argument delimit for each file  
#### `usehead`:: The header from which input file will be used for output file.  
	- Default: None (Don't write header)
#### `gzip`:: argument gzip for each file  
#### `match`:: The match function.   
#### `do`:: The do function. Global vaiable `fout` is available to write results to output file.  

### requires
[`python-simread`](https://github.com/pwwang/simread)
{% endraw %}
