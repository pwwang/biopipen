# eqtl
<!-- toc -->
{% raw %}

## pMatrixeQTL

### description
Call eQTLs using Matrix eQTL

### input
#### `snpfile:file`:: The genotype file, rows are snps and columns are samples  
#### `expfile:file`:: The expression file, rows are genes  
#### `covfile:file`:: The covariant file, rows are covariants  

### output
#### `outfile:file`:: The matrix eqtl output file  

### args
#### `model`:: The model to use, either modelLINEAR(default) or modelANOVA  
#### `pval` :: The pvalue cutoff (if `cisopts.dist` > 0, will be used as pval for trans-eQTL)  
#### `fdr`  :: Calculate FDR or not (default: True)  
#### `cisopts`:: Options for calling cis-, trans-eQTL  
	- `snppos` : The snp position file (columns are: snp, chr, pos)
	- `genepos`: The gene position file (columns are: gene, chr, start, end)
	- `dist`   : The distance to define cis-eQTL. (default: 0 (don't do cis-, trans- calling)
	- `cispv`  : The pvalue cutoff for cis-eQTL (`pval` will not work)
{% endraw %}
