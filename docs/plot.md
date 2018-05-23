# plot
<!-- toc -->
{% raw %}

## pPlot

### description
Use ggplot2 to generate plots

### infile
#### `infile:file`:: The input data file  

### outfile
#### `outfile:file`:: The output file  

### args
#### `cnames` :: Whether the input file has colnames. Default: True  
#### `rnames` :: Whether the input file has rownames. Default: False  
#### `aes`    :: The default aes. Default: {'x':1, 'y':2} (corresponding to colnames)  
#### `helper` :: Some helper codes to generate `params` and `ggs`  
#### `devpars`:: The device parameters. Default: `Box(res = 300, height = 2000, width = 2000)`  
#### `ggs`    :: The extra ggplot elements.  

## pScatter

### description
Use ggplot2 geom_point to generate plots

### infile
#### `infile:file`:: The input data file  

### outfile
#### `outfile:file`:: The output file  

### args
#### `cnames` :: Whether the input file has colnames. Default: True  
#### `rnames` :: Whether the input file has rownames. Default: False  
#### `x`      :: The x aes. Default: 1 (corresponding to colnames)  
#### `y`      :: The y aes. Default: 2 (corresponding to colnames)  
#### `helper` :: Some helper codes to generate `params` and `ggs`  
#### `devpars`:: The device parameters. Default: `Box(res = 300, height = 2000, width = 2000)`  
#### `params` :: The extra params for `geom_point`  
#### `ggs`    :: The extra ggplot elements.  

## pPoints

### description
Alias for pScatter

## pHisto

### description
Use ggplot2 geom_histogram to generate histograms

### infile
#### `infile:file`:: The input data file  

### outfile
#### `outfile:file`:: The output file  

### args
#### `cnames` :: Whether the input file has colnames. Default: True  
#### `rnames` :: Whether the input file has rownames. Default: False  
#### `x`      :: The x aes. Default: 1 (corresponding to colnames)  
#### `helper` :: Some helper codes to generate `params` and `ggs`  
#### `devpars`:: The device parameters. Default: `Box(res = 300, height = 2000, width = 2000)`  
#### `params` :: The extra params for `geom_point`  
#### `ggs`    :: The extra ggplot elements.  

## pFreqpoly

### description
Use ggplot2 geom_freqpoly to generate frequency polygon plot.

### infile
#### `infile:file`:: The input data file  

### outfile
#### `outfile:file`:: The output file  

### args
#### `cnames` :: Whether the input file has colnames. Default: True  
#### `rnames` :: Whether the input file has rownames. Default: False  
#### `x`      :: The x aes. Default: 1 (corresponding to colnames)  
#### `helper` :: Some helper codes to generate `params` and `ggs`  
#### `devpars`:: The device parameters. Default: `Box(res = 300, height = 2000, width = 2000)`  
#### `params` :: The extra params for `geom_point`  
#### `ggs`    :: The extra ggplot elements.  

## pBoxplot

### description
Generate box plot

### input
#### `datafile:file`:: The data file  

### output
#### `outpng:file`:: The output figure  

### args
#### `header`::    `header` parameter for `read.table`, default: True  
#### `rownames`::  `row.names` parameter for `read.table`, default: 1  
#### `params`::    Other parameters for `boxplot`, default: ""  

## pHeatmap

### description
Plot heatmaps.

### input
#### `infile:file`:: The input matrix file  

### output
#### `outfile:file`:: The heatmap  

### args
#### `ggs`:: The ggplot items for heatmap  
#### `devpars`:: The parameters for device. Default: `{'res': 300, 'height': 2000, 'width': 2000}`  
#### `dendro`:: The parameters for control of the dendrogram. Default: `{'dendro': True}`  
	- `dendro`: `True`: plot dendros for both rows and cols; `col`: only plot dendro for cols; `row`: only plot dendro for rows
	- `rows`: The rownames to subset the rows and control the order of rows. Must a list. Only works when not plotting dendro for rows.
	- `cols`: The colnames to subset the cols and control the order of cols. Must a list. Only works when not plotting dendro for cols.
#### `header`:: The input file has header? Default: True  
#### `rownames`:: The input file has rownames? Default: 1  
#### `rows`:: Row selector  
	- `all`: All rows
	- `top:N`: Top N rows (original data ordered in descending order). N defaults to 100
	- `bottom:N`: Bottom N rows. N defaults to 100
	- `both:N`: Top N rows and bottom N rows. N defaults to 50
	- `random:N`: Random N rows. N defaults to 50
	- `random-both:N`: Random N rows from top part and N rows from bottom part. N defaults to 50
#### `cols`:: Col selector (see `rows`).  

## pScatterCompare

### description
Plot scatter plot to compare values of first 2 columns of input data

### input
#### `infile:file`:: The input file containing a matrix with at least 2 columns  
	- Other columns are groups used to group the scatter points
	- Data must be normalized to [0, 1]

### output
#### `outfile:file`:: The output plot  

### args
#### `ggs`:: Extra expressions for ggplot. Note if geom_point is included, original geom_point will be ignored.  
#### `devpars`:: The parameters for plot device. Default: `{'res': 300, 'height': 2000, 'width': 2000}`  
#### `rownames`:: Whether the input file has row names. Default: True  
#### `regr`:: Whether draw the regression line. Default: False  
#### `corr`:: The method to calculate the correlation. Default: `pearson`  
	- Could be: `pearson`, `spearman` or `kendall`
	- If it's neither of the three, no correlations will show.

## pROC

### description
Generate ROC curves and output AUC.

### input
#### `infile:file`:: The input matrix file.  
	- Col1: rownames if args.rnames is True else label (0, 1 class)
	- Col2: prediction values from model1
	- ...

## pVenn

### description
Venn/UpsetR plots.

### input
#### `infile:file`:: The input matrix  
	- format:
	```
		category1	category2	category3
	[e1]	0	1	1
	[e2]	0	0	1
	...
	[eN]	1	0	0
	```
	rownames are not necessary but colnames are.

### output
#### `outfile:file`:: The plot  

### args
#### `tool`    :: Which tools to use. Default: auto (venn, upsetr, auto(n<=3: venn, otherwise upsetr))  
#### `rnames`  :: Whether input file has rownames. Default: False  
#### `params`  :: Other params for `venn.diagram` or `upset`. Default: {}  
#### `devpars` :: The parameters for plot device. Default: `{'res': 300, 'height': 2000, 'width': 2000}`  

### requires
[`r-VennDiagram`](https://www.rdocumentation.org/packages/VennDiagram)
[`r-UpSetR`](https://www.rdocumentation.org/packages/UpSetR)

## pPie

### description
Plot piechart

### input
#### `infile:file`:: The input file. Could be either:  
	- Direct numbers of each category.
	```
	Group1	Group2
	50	50
	```
	- Presence of each items in the category.
	```
		Group1	Group2
	Item1	1	0
	Item2	0	1
	...
	```

### output
#### `outfile:file`:: the output plot  

### args
#### `rnames` :: Whether the input file has row names. Default: `False`  
#### `ggs`    :: Extra expressions for ggplot.  
#### `devpars`:: The parameters for plot device. Default: `{'res': 300, 'height': 2000, 'width': 2000}`  
{% endraw %}
