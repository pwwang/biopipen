# rank
<!-- toc -->
{% raw %}

## pRankProduct

### description
Calculate the rank product of a set of ranks. Refer to [here](https://en.wikipedia.org/wiki/Rank_product)

### input
#### `infile:file`:: The input file  
- Format:
```
			Case1	Case2	...
Feature1	8.2  	10.1 	...
Feature2	2.3  	8.0  	...
...
```
- Or instead of values, you can also have ranks in the input file:
```
			Rank1	Rank2	...
Feature1	2    	1    	...
Feature2	3    	2    	...
...
```

### output
#### `outfile:file`:: The output file with original ranks, rank products and p-value if required  

### args
#### `informat`:: The input format of the values. Whether they are real values (value) or ranks (rank). Default: value  
#### `pval`::     Whether to calculate the p-value or not. Default: True  
#### `header`::   Whether the input file has headers (rownames are required!). Default: True  
#### `plot`::     Number of rows to plot. Default: 0 (Don't plot)  
#### `cex`::      Font size for plotting. Default: 0.9  
#### `cnheight`:: Colname height. Default: 80  
#### `rnwidth`::  Rowname width. Default: 50  
#### `width`::    Width of the png file. Default: 2000  
#### `height`::   height of the png file. Default: 2000  
{% endraw %}
