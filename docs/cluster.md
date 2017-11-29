# cluster
<!-- toc -->
{% raw %}

## pDist2Coords

### description
Convert a distance matrix to coordinates, using multidimensional scaling.

### input
#### `infile:file`:
The distance matrix, could be a full distance matrix, a triangle matrix or a pair-wise distance file  
	- full dist matrix (full):
	```
		s1	s2	s3
	s1	0	1	1
	s2	1	0	1
	s3	1	1	0
	```
	- triangle matrix (upper/lower), could be also lower triangle
	```
		s1	s2	s3
	s1	0	1	1
	s2		0	1
	s3			0
	```
	- pair-wise (pair): (assuming auto-pair-wise distance = 0, that is: `s1	s1	0`)
	```
	s1	s2	1
	s1	s3	1
	s2	s3	1
	```
	- Both rownames and header are required.

### output
#### `outfile:file`:
The output coordinate file  

## pCluster

### description
Use `optCluster` to select the best number of clusters and cluster method, then perform the clustering

### input
#### `infile:file`:
The input matrix file. Clustering will be performed against rows. If not, set `args.transpose` = True  

### output
#### `outfile:file`:
The output cluster file  
#### `outdir:dir`:
The output directory containing the figures  

### args
#### `transpose`:
Transpose the input matrix. Default: False  
#### `cnames`:
Whether the input matrix contains header before transposing. Default: False  
#### `rnames`:
Which column is the rownames before transposing. Default: 1  
#### `plot`:
Whether plot the cluster. Default: True  
#### `minc`:
Min number of clusters to test. Default: 2  
#### `maxc`:
Min number of clusters to test. Default: 15  
	- If number of rows (nrows) <= 15, then max = nrows - 1
#### `methods`:
The methods to test. Default: "all"  
	- Could be any of "agnes", "clara", "diana", "fanny", "hierarchical", "kmeans", "model", "pam", "som", "sota", "em.nbinom", "da.nbinom", "sa.nbinom", "em.poisson", "da.poisson", "sa.poisson"
	- Multiple methods could be separated by comma (,), or put in a list
	- By default, fanny, model and sota will be excluded because fanny causes error and the latter two are slow. You can manually include them if you want.
	- Improper methods will be automatically excluded by `args.isCount`
#### `isCount`:
Whether the data is count data. Corresponding methods will be tested. Default: False  

## pMCluster

### description
Use `r-mclust` to do clustering. Current just do simple clustering with the package

### input
#### `infile:file`:
The input a coordinate file  

### output
#### `outdir:dir`:
The output of final results  

### args
#### `transpose`:
Transpose the input matrix? Default: False  
#### `rnames`:
The `row.names` for `read.table` to read the input file, default: True.  
#### `cnames`:
The `header` argument for `read.table` to read the input file, default: True.  
#### `caption`:
The caption for the `fviz_cluster`, default: "CLARA Clustering".  
#### `minc`:
The min # clusters to try, default: 2  
#### `maxc`:
The max # clusters to try, default: 15  

## pAPCluster

### description
Use `r-apcluster` to do clustering. 

### input
#### `infile:file`:
The input a coordinate file  

### output
#### `outdir:dir`:
The output of final results  

### args
#### `transpose`:
Transpose the input matrix? Default: False  
#### `rownames`:
The `row.names` for `read.table` to read the input file, default: 1.  
#### `header`:
The `header` argument for `read.table` to read the input file, default: True.  
#### `caption`:
The caption for the `fviz_cluster`, default: "APClustering".  

## pHCluster

### description
Do hierarchical clustering.

### input
#### `infile:file`:
The input files with variants as rows, features as columns.  
	- NOTE: clustering is performed on rows, rownames are the leaf labels.

### output
#### `outdir:dir`:
The result directory, containing:  
	- `hclust.merge.txt`: including merge and height information
	- `hclust.order.txt`: including order and labels information
	- `hclust.png`:       the dendrogram plot

### args
#### `fast`:
whether to use `fastcluster` package or not, default: False  
#### `gg`:
whether to use `ggdendro` or not, default: False  
#### `rownames`:
The `row.names` for `read.table` to read the input file, default: 1.  
#### `header`:
The `header` argument for `read.table` to read the input file, default: True.  
#### `method`:
Which method to use for `hclust`. Default: "complete" (use `?hclust` to check all availables)  
#### `rotate`:
Which to rotate the plot or not. Default: False  
#### `transpose`:
Whether to transpose the matrix before cluster. Default: False  
{% endraw %}
