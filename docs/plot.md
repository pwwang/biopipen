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

## pScatter

### description
Use ggplot2 geom_point to generate plots

### infile
#### `infile:file`:: The input data file  

### outfile
#### `outfile:file`:: The output file  

## pPoints

## pHisto

### description
Use ggplot2 geom_histogram to generate histograms

### infile
#### `infile:file`:: The input data file  

### outfile
#### `outfile:file`:: The output file  

## pFreqpoly

### description
Use ggplot2 geom_freqpoly to generate frequency polygon plot.

### infile
#### `infile:file`:: The input data file  

### outfile
#### `outfile:file`:: The output file  

## pBoxplot

### description
Generate box plot

### input
#### `datafile:file`:: The data file  

### output
#### `outpng:file`:: The output figure  

## pHeatmap

### description
Plot heatmaps.

### input
#### `infile:file`:: The input matrix file  

### output
#### `outfile:file`:: The heatmap  

## pScatterCompare

### description
Plot scatter plot to compare values of first 2 columns of input data

### input
#### `infile:file`:: The input file containing a matrix with at least 2 columns  
	- Other columns are groups used to group the scatter points
	- Data must be normalized to [0, 1]

### output
#### `outfile:file`:: The output plot  

## pROC

### description
Generate ROC curves and output AUC.

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
{% endraw %}
