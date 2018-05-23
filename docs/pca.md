# pca
<!-- toc -->
{% raw %}

## pPCA

### description
Perform PCA analysis

### input
#### `infile:file`:: The matrix to do the analysis  
- Note that rows are samples, columns are features, if not, use `args.transpose = True`

### output
#### `outfile:file`:: The output coordinate file  
- Columns are PCs, rows are samples

### args
#### `transpose`::  Whether to transpose the input matrix from infile. Default: False  
#### `rownames`::   The `row.names` argument for `read.table`, default: 1  
#### `header`::     The `header` argument for `read.table` to read the input file, default: True.  
#### `screeplot`::  Whether to generate the screeplot or not. Default: True  
#### `sp_ncp`::     Number of components in screeplot. Default: 0 (auto detect)  
- if total # components (tcp) < 20: use all
- else if tcp > 20, use 20
#### `varplot`::    Whether to generate the variable plot or not. Default: False  
#### `biplot`::     Whether to generate the variable plot or not. Default: True  

### requires
[`r-factoextra`](https://cran.r-project.org/web/packages/factoextra/index.html) for plots

## pSelectPCs

### description
Select a subset of PCs from pPCA results

### input
#### `indir:file`:: The directory generated from pPCA  

### output
#### `outfile:file`:: The file containing selected PCs  

### args
#### `n`:: The number of PCs to select. Default: 0.9  
- If it is < 1, used as the % variation explained from stdev.txt
{% endraw %}
