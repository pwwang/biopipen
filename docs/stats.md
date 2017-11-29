<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**  *generated with [DocToc](https://github.com/thlorenz/doctoc)*

- [pMetaPval](#pmetapval)
- [pMetaPval1](#pmetapval1)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->


## pMetaPval

### description
	Combine p-values in the files from input directory

### input
#### `indir:dir`:
 The directory containing the input files  

### output
#### `outfile:file`:
 The output file containing the meta-pvalues  

### args
#### `args.pattern`:
 The pattern used to filter the input files. Default: '*'  
#### `args.header`:
 Whether the input files contains a header. Default: True  
		- Could be a list to specify it for each file.
		- The order should be concordant with the file names
#### `args.pcol`:
 Which column is the p-value. Default: -1 (last column)  
#### `args.poutonly`:
 Only output pvalues. Default: False (output all possible information)  
#### `args.outheader`:
 Whether output the header. Default: True  
#### `args.method`:
 The method used to calculate the meta-pvalue. Default: sumlog (Fisher's method)  
		- Other available methods: logitp, sumz, votep, sump, meanp and wilkinsonp
		- See: https://www.rdocumentation.org/packages/metap/versions/0.8

## pMetaPval1

### description
	Combine p-values in a single file by rows.

### input
#### `infile:file`:
 The input file  

### output
#### `outfile:file`:
 The output file containing the meta-pvalues  

### args
#### `args.header`:
 Whether the input files contains a header. Default: True  
#### `args.pcol`:
 Which column is the p-value. Default: -1 (last column)  
#### `args.poutonly`:
 Only output pvalues. Default: False (output all possible information)  
#### `args.outheader`:
 Whether output the header. Default: True  
#### `args.method`:
 The method used to calculate the meta-pvalue. Default: sumlog (Fisher's method)  
		- Other available methods: logitp, sumz, votep, sump, meanp and wilkinsonp
		- See: https://www.rdocumentation.org/packages/metap/versions/0.8
