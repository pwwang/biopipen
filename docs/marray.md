# marray
<!-- toc -->
{% raw %}

## pCELdir2Matrix

### description
Convert CEL files to expression matrix
File names will be used as sample names (colnames)

### input
#### `indir:file`::  the directory containing the CEL files, could be gzipped  
	- If you have files, then use `pFiles2Dir` first

### output
#### `outfile:file`:: the expression matrix file  
#### `outdir:dir`::   the directory containing expr file and plots  

### args
#### `pattern`  :: The pattern to filter files. Default `'*'`  
#### `norm`     :: The normalization method. Default: rma (mas5)  
#### `gfile`    :: The group file. Default: ''  
#### `cdffile`  :: The cdffile. Default: ''  
#### `annofile` :: The annotation file. Default: ''  
#### `hmrows`   :: How many rows to be used to plot heatmap  
#### `plot`:: Whether to plot  
	- `boxplot`   : Whether to plot a boxplot. Default: False
	- `heatmap`   : Whether to plot a heatmap. Default: False
	- `histogram` : Whether to plot a histgram. Default: False
#### `devpars`    :: Parameters for png. Default: `{'res': 300, 'width': 2000, 'height': 2000}`  
#### `ggs`:: The ggplot parameters  
	- `boxplot`  : The ggplot parameters for boxplot. Default: `Box(ylab = {0: "Log2 Intensity"})`
	- `heatmap`  : The ggplot parameters for heatmap. Default: `Box(theme = {'axis.text.y': 'r:element_blank()'})`
	- `histogram`: The ggplot parameters for histgram. Default: `Box(labs = {'x': "Log2 Intensity", "y": "Density"})`
{% endraw %}
