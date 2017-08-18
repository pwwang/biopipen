from pyppl import proc

"""
@name:
	pDist2Coords
@description:
	Convert a distance matrix to 2D coordinates, using multidimensional scaling
@input:
	`infile:file`: The distance matrix, could be a full distance matrix, a triangle matrix or a pair-wise distance file
	- full dist matrix (full):
	```
		s1	s2	s3
	s1	0	1	1
	s2	1	0	1
	s3	1	1	0
	```
	- triangle matrix (triangle), could be also lower triangle
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
	- Both rownames and header of `full` and `triangle` can be omitted, just set `args.rownames = "NULL"` and `args.header = False`
@output:
	`outfile:file`: The output coordinate file
@args:
	`informat`: The format of the input file: full, triangle or pair. Default: full
	`rownames`: The `row.names` argument for `read.table`, default: 1
	`header`:   The `header` argument for `read.table` to read the input file, default: True.
	`k`:        How many dimension? Default: 2 (R^2)
"""
pDist2Coords = proc()
pDist2Coords.input   = "infile:file"
pDist2Coords.output  = "outfile:file:{{infile | fn}}.coords"
pDist2Coords.args    = {'informat': 'full', 'rownames': 1, 'header': True, 'k': 2}
pDist2Coords.lang    = "Rscript"
pDist2Coords.script  = """
informat = "{{args.informat}}"
rownames = {{args.rownames}}
if (informat == "full") {
	mat = read.table ("{{infile}}", row.names=rownames, header={{args.header | Rbool }}, check.names=F)
} else if (informat == "triangle") {
	mat = read.table ("{{infile}}", sep="\\t", row.names=rownames, header={{args.header | Rbool }}, check.names=F, fill=T)
	for (i in 1:nrow(mat)) {
		for (j in i:nrow(mat)) {
			if (j == i)  { mat[i,j] = 0; mat[j, i] = 0}
			else if (is.na(mat[i,j])) mat[i,j] = mat[j,i]
			else mat[j,i] = mat[i,j]
		}
	}
} else {
	f <- read.table("{{infile}}", header=F, row.names=NULL, stringsAsFactor=F, check.names=F)
	f[, 1] = as.character (f[, 1])
	f[, 2] = as.character (f[, 2])
	name1 = f[,1]
	name2 = f[,2]
	names = union(name1, name2)
	nlen  = length(names)
	mat   = matrix(NA, ncol=nlen, nrow=nlen)
	colnames(mat) = names
	rownames(mat) = names
	for (i in 1:nrow(f)) {
		mat[f[i,1], f[i,1]] = 0
		mat[f[i,2], f[i,2]] = 0
		mat[f[i,1], f[i,2]] = f[i,3]
		mat[f[i,2], f[i,1]] = f[i,3]
	}
	print (mat)
}
mat[is.na(mat)] = sum(mat)
coords = cmdscale(mat, k={{args.k}})
write.table (coords, "{{outfile}}", col.names=F, row.names=F, quote=F, sep="\\t")
"""

"""
@name:
	pCluster
@description:
	Use `optCluster` to select the best number of clusters and cluster method, then perform the clustering
@input:
	`infile:file`: The input matrix file. Clustering will be performed against rows. If not, use `args.transpose` = True
@output:
	`outfile:file`: The output cluster file
	`outdir:dir`:   The output directory containing the figures
@args:
	`transpose`:    Transpose the input matrix. Default: False
	`header`:       Whether the input matrix contains header before transposing. Default: False
	`rownames`:     Which column is the rownames before transposing. Default: 1
	`plot`:         Whether plot the cluster. Default: True
	`nc`:           Number of clusters to test. Default: "2:15"
	`methods`:      The methods to test. Default: "all"
		- Could be any of "agnes", "clara", "diana", "fanny", "hierarchical", "kmeans", "model", "pam", "som", "sota", "em.nbinom", "da.nbinom", "sa.nbinom", "em.poisson", "da.poisson", "sa.poisson"
		- Multiple methods could be separated by comma (,), or put in a list
		- By default, fanny, model and sota will be excluded because fanny causes error and the latter two are slow. You can manually include them if you want.
		- Improper methods will be automatically excluded by `args.isCount`
	`isCount`:      Whether the data is count data. Corresponding methods will be tested. Default: False
@requires:
	[`r-optCluster`](https://rdrr.io/cran/optCluster/man/optCluster.html)
	[`r-factoextra`](https://cran.r-project.org/web/packages/factoextra/index.html)
"""
pCluster                = proc (desc = 'Use `optCluster` to select the best number of clusters and cluster method, then perform the clustering.')
pCluster.input          = 'infile:file'
pCluster.output         = 'outfile:file:{{infile | fn}}.cluster/clusters.txt, outdir:dir:{{infile | fn}}.cluster'
pCluster.args.transpose = False
pCluster.args.header    = False
pCluster.args.rownames  = 1
pCluster.args.plot      = True
pCluster.args.nc        = "2:15"
# "agnes", "clara", "diana", "fanny", "hierarchical", "kmeans", "model", "pam", "som", "sota", "em.nbinom", "da.nbinom", "sa.nbinom", "em.poisson", "da.poisson", "sa.poisson"
# Note fanny reports error and the algorithm cannot go further.
# It is excluded by default, you can add it manually
# model and sota are excluded because they are slow
# You can also manually add them
pCluster.args.methods   = 'all'
pCluster.args.isCount   = False
pCluster.lang           = 'Rscript'
pCluster.script         = """
library(optCluster)

d = read.table ("{{infile}}", sep="\\t", header={{args.header | Rbool}}, row.names={{args.rownames}}, check.names=F)
if ({{args.transpose | Rbool}}) {
	d = t(d)
}
names = rownames(d)
methods = {{args.methods | lambda x: 'c("'+ '","'.join(x) +'")' if isinstance(x, list) else ('"all"' if x == 'all' else 'c("'+ '","'.join(map(lambda y: y.strip(), x.split(','))) +'")')}}
if (methods == 'all') {
	methods = c("agnes", "clara", "diana", "hierarchical", "kmeans", "pam", "som", "em.nbinom", "da.nbinom", "sa.nbinom", "em.poisson", "da.poisson", "sa.poisson")
}
if ({{args.isCount | Rbool}}) {
	methods = intersect(c("em.nbinom", "da.nbinom", "sa.nbinom", "em.poisson", "da.poisson", "sa.poisson"), methods)
} else {
	methods = intersect(c("agnes", "clara", "diana", "hierarchical", "kmeans", "model", "pam", "som", "sota", "fanny"), methods)
}
ret = optCluster(d, {{args.nc}}, clMethods = methods, countData = {{args.isCount | Rbool}}, clVerbose = T, rankVerbose = T)
ret = optAssign(ret)

outfile = "{{outdir}}/clusters.txt"
metfile = "{{outdir}}/method.txt"
outfig  = "{{outdir}}/clusters.png"
outdir  = "{{outdir}}"
clusters = matrix(ret$cluster, ncol=1)
rownames(clusters) = names
write.table (clusters, outfile, row.names=T, sep="\\t", quote=F, col.names=F)
write.table (ret$method, metfile, row.names=F, sep="\\t", quote=F, col.names=F)
if ({{args.plot | Rbool}}) {
	library(factoextra)
	pdata = list()
	pdata$data = d
	pdata$cluster = ret$cluster
	png (file = outfig, res=300, width=2000, height=2000)
	print (fviz_cluster(pdata, data=d))
	dev.off()
}
"""

"""
@name:
	pMClust
@description:
	Use `r-mclust` to do clustering. Current just do simple clustering with the package
@input:
	`infile:file`:  The input a coordinate file
@output:
	`outdir:dir`:   The output of final results
@args:
	`rownames`:     The `row.names` for `read.table` to read the input file, default: 1.
	`header`:       The `header` argument for `read.table` to read the input file, default: True.
	`caption`:      The caption for the `fviz_cluster`, default: "CLARA Clustering".
	`min`:          The min # clusters to try, default: 2
	`max`:          The max # clusters to try, default: 15
@requires:
	[`r-mclust`](https://cran.r-project.org/web/packages/mclust/index.html)
	[`r-factoextra`](https://cran.r-project.org/web/packages/factoextra/index.html)
"""
pMClust  = proc()
pMClust.input   = "infile:file"
pMClust.output  = 'outdir:dir:{{infile | fn}}.mclust'
pMClust.args    = { 'rownames': 1, 'header': True, 'min':2, 'max':15, 'caption': 'MClust Clustering' }
pMClust.lang    = 'Rscript'
pMClust.script  = """
library ('factoextra')
library ('mclust')
data = read.table ("{{infile}}", row.names={{args.rownames}}, header={{args.header | Rbool }}, check.names=F)
data = data[, apply(data, 2, sd)!=0, drop=F]
mc        = Mclust (data, G={{args.min}}:{{args.max}})
clustfile = file.path ("{{outdir}}", "cluster.txt")
write.table (mc$classification, clustfile, quote=F, col.names = F, sep="\\t")
clustfig  = file.path ("{{outdir}}", "cluster.png")
png (file = clustfig)
labels    = colnames(data)
print (fviz_cluster (mc, main = gsub("%K%", max(mc$classification), "{{args.caption}}")))
dev.off()
"""

"""
@name:
	pAPCluster
@description:
	Use `r-apcluster` to do clustering. 
@input:
	`infile:file`:  The input a coordinate file
@output:
	`outdir:dir`:   The output of final results
@args:
	`rownames`:     The `row.names` for `read.table` to read the input file, default: 1.
	`header`:       The `header` argument for `read.table` to read the input file, default: True.
	`caption`:      The caption for the `fviz_cluster`, default: "APClustering".
@requires:
	[`r-apcluster`](https://cran.r-project.org/web/packages/apcluster/index.html)
	[`r-factoextra`](https://cran.r-project.org/web/packages/factoextra/index.html)
"""
pAPCluster  = proc()
pAPCluster.input   = "infile:file"
pAPCluster.output  = 'outdir:dir:{{infile | fn}}.apcluster'
pAPCluster.args    = { 'rownames': 1, 'header': True, 'caption': 'APClustering with K=%K%' }
pAPCluster.lang    = 'Rscript'
pAPCluster.script  = """
library ('factoextra')
library ('apcluster')
data = read.table ("{{infile}}", row.names={{args.rownames}}, header={{args.header | Rbool }}, check.names=F)
data = data[, apply(data, 2, sd)!=0, drop=F]
ap        = apcluster(negDistMat(r=2), data)
clusters  = ap@clusters
nclust    = length (clusters)
clusts    = matrix (NA, ncol=1, nrow=ap@l)
rownames (clusts) = rownames(data)
for (i in 1:nclust) {
	cids  = unlist (clusters[i])
	for (cid in cids) clusts[cid, 1] = i
}
fvizobj   = {}
fvizobj$data = data
fvizobj$cluster = clusts

clustfile = file.path ("{{outdir}}", "cluster.txt")
write.table (clusts, clustfile, quote=F, col.names = F, sep="\\t")
clustfig  = file.path ("{{outdir}}", "cluster.png")
png (file = clustfig)
labels    = colnames(data)
print (fviz_cluster (fvizobj, main = gsub("%K%", nclust, "{{args.caption}}")))
dev.off()
"""

"""
@name:
	pHClust
@description:
	Do hierarchical clustering.
@input:
	`infile:file`: The input files with variants as rows, features as columns.
	- NOTE: clustering is performed on rows, rownames are the leaf labels.
@output:
	`outdir:dir`:  The result directory, containing:
	- `hclust.merge.txt`: including merge and height information
	- `hclust.order.txt`: including order and labels information
	- `hclust.png`:       the dendrogram plot
@args:
	`fast`:     whether to use `fastcluster` package or not, default: False
	`gg`:       whether to use `ggdendro` or not, default: False
	`rownames`: The `row.names` for `read.table` to read the input file, default: 1.
	`header`:   The `header` argument for `read.table` to read the input file, default: True.
	`method`:   Which method to use for `hclust`. Default: "complete" (use `?hclust` to check all availables)
	`rotate`:   Which to rotate the plot or not. Default: False
	`transpose`:Whether to transpose the matrix before cluster. Default: False
@requires:
	[`r-fastcluster`](https://cran.r-project.org/web/packages/fastcluster/index.html) if `args.fast` is True
	[`r-ggdendro`](https://cran.r-project.org/web/packages/ggdendro/index.html) if `args.gg` is True
"""
pHClust = proc()
pHClust.input  = "infile:file"
pHClust.output = "outdir:dir:{{infile | fn}}.hclust"
pHClust.args   = {"fast":False, "gg":False, "rownames":1, "header":True, 'method': 'complete', 'rotate': False, 'transpose': False}
pHClust.lang   = "Rscript"
pHClust.script = """
data = read.table ("{{infile}}", row.names={{args.rownames}}, header={{args.header | Rbool }}, check.names=F)
if ({{args.transpose | Rbool}}) data = t(data)
data = data[, apply(data, 2, sd)!=0, drop=F]
dmat = dist(data)
if ({{args.fast | Rbool}}) {
	library('fastcluster')
}
orderfile = file.path ("{{outdir}}", "hclust.order.txt")
mergefile = file.path ("{{outdir}}", "hclust.merge.txt")
clustfig  = file.path ("{{outdir}}", "hclust.png")
hobj = hclust (dmat, method="{{args.method}}")
png (file = clustfig)
if ({{args.gg | Rbool}}) {
	library('ggplot2')
	library('ggdendro')
	ggdendrogram (hobj, rotate={{args.rotate | Rbool}})
} else {
	plot (as.dendrogram(hobj), horiz={{args.rotate | Rbool}})
}
dev.off()
orderobj = data.frame (Order=hobj$order, Labels=hobj$labels)
mergeobj = data.frame (Merge1=hobj$merge[,1], Merge2=hobj$merge[,2], Height=hobj$height)
write.table (orderobj, orderfile, quote=F, row.names=F, col.names = T, sep="\\t")
write.table (mergeobj, mergefile, quote=F, row.names=F, col.names = T, sep="\\t")

"""
