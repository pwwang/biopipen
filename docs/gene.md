# gene
<!-- toc -->
{% raw %}

## pGeneNameNorm

### description
Normalize gene names using MyGeneinfo.

### input
#### `infile:file`:: The input file  

### output
#### `outfile:file`:: The output file  

### args
#### `notfound`:: What if a symbol is not found. Default: ignore  
	- skip  : skip the record(don't write it to output file)
	- ignore: use the original name;
	- error : report error
#### `col`:: the column index containing the gene names  
#### `from`:: the original format. Default: 'symbol, alias'  
#### `to`:: the output gene name format. Default: 'symbol'  
#### `genome`:: the genome. Default: 'hg19'  

## pGeneTss

### description
Get gene TSS in BEd format.

### input
#### `infile:file`:: The input file containing genes  

### output
#### `outfile:file`:: The output BED file  

### args
#### `notfound`:: What if the gene is not found. Default: skip.  
	- error: report error
#### `header`:: Whether the input file contains header. Default: False  
#### `skip`:: Skip N lines of input file. Default: 0  
	- This has highest priority of header and comment
#### `comment`:: The comment line start sign. Default: #  
#### `delimit`:: The delimit of input file if it has multiple column. Default: `\\t`  
#### `col`:: The column index contains the genes. Default: 0  
#### `frm`:: The format of the genes. Default: `symbol, alias`  
#### `tmpdir`:: The tmpdir used to store mygene cache files.  
#### `genome`:: In which genome to fetch the coordinates. Default: hg19  

## pGeneBody

### description
Get gene body region in BED format

### input
#### `infile:file`:: The input file containing genes  

### output
#### `outfile:file`:: The gene body region  

### args
#### `notfound`:: What if a gene is not found when transfer the gene names to gene symbols  
	- error: report error
	- skip (default): skip it
#### `inmeta`::   The metadata for input file, mainly to indicate where the GENE column is.  
#### `inopts`::   Input options for reading input file.  
	- skip: number of lines to skip. Default: 0
	- comment: the starting string for comment lines. Default: #
	- delimit: The delimit for the input file. Default: '\\t'
frm: The gene name format in the input file. Default: 'symbol, alias'
tmpdir: The tmpdir to cache the gene name conversion.
genome: The genome used to do the conversion.
{% endraw %}
