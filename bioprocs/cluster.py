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
	- Both rownames and header of `full` and `triangle` can be omitted, just set `proc.args.rownames = "NULL"` and `proc.args.header = False`
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
informat = "{{proc.args.informat}}"
rownames = {{proc.args.rownames}}
if (informat == "full") {
	mat = read.table ("{{infile}}", row.names=rownames, header={{proc.args.header | Rbool }}, check.names=F)
} else if (informat == "triangle") {
	mat = read.table ("{{infile}}", sep="\\t", row.names=rownames, header={{proc.args.header | Rbool }}, check.names=F, fill=T)
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
coords = cmdscale(mat, k={{proc.args.k}})
write.table (coords, "{{outfile}}", col.names=F, row.names=F, quote=F, sep="\\t")
"""

"""
@name:
	pDecideK
@description:
	Decide number of clusters using different methods
@input:
	`infile:file`: the coordinates file, if all you have is a distance/similarity file, convert it to coordinates file using `pDist2Coords`
@output:
	`kfile:file`: the output file with `K`
@args:
	`method`:                         The method used to determine K
	- `elbow:<ev.thres>:<inc.thres>`: Look for a bend or elbow in the sum of squared error (SSE) scree plot, see [ref](https://artax.karlin.mff.cuni.cz/r-help/library/GMD/html/elbow.html). Default: `elbow` = `elbow:.95:01`
	- `pamk:<min>:<max>`:             You can do partitioning around medoids to estimate the number of clusters using the pamk function in the fpc package. Default: `pamk` = `pamk:2:15`
	- `calinski:<min>:<max>`:         Calinski criterion. Default: `calinski` means `calinski:2:15`
	- `mclust:<min>:<max>`:           Determine the optimal model and number of clusters according to the Bayesian Information Criterion for expectation-maximization, initialized by hierarchical clustering for parameterized Gaussian mixture models. [Ref1](http://www.stat.washington.edu/research/reports/2006/tr504.pdf
#), [Ref2](http://www.jstatsoft.org/v18/i06/paper). Default: `mclust` = `mclust:2:15`
	- `ap`:                           Affinity propagation (AP) clustering, see [ref](http://dx.doi.org/10.1126/science.1136800)
	- `gap:<min>:<max>`:              Gap Statistic for Estimating the Number of Clusters. Default: `gap:2:10`
	- `nbclust`:                      The [NbClust package](http://cran.r-project.org/web/packages/NbClust/index.html) provides 30 indices to determine the number of clusters in a dataset.
	`rownames`:                       The `row.names` for `read.table` to read the input file, default: 1.
	`header`:                         The `header` argument for `read.table` to read the input file, default: True.
	`seed`:                           The seed for randomization, default: 0.
@requires:
	[`r-cluster`](https://cran.r-project.org/web/packages/cluster/index.html) if `gap` method used
	[`r-GMD`](https://cran.r-project.org/web/packages/GMD/index.html) if `elbow` method userd
	[`r-fpc`](https://cran.r-project.org/web/packages/fpc/index.html) if `pamk` method used
	[`r-vegan`](https://cran.r-project.org/web/packages/vegan/index.html) if `calinski` method used
	[`r-mclust`](https://cran.r-project.org/web/packages/mclust/index.html) if `mclust` method used
	[`r-apcluster`](https://cran.r-project.org/web/packages/apcluster/index.html) if `ap` method used
	[`r-NbClust`](https://cran.r-project.org/web/packages/NbClust/index.html) if `nbclust` method used
"""
pDecideK = proc()
pDecideK.input   = "infile:file"
pDecideK.output  = "kfile:file:{{infile | fn}}-K.txt"
pDecideK.args    = { 'method': 'elbow', 'rownames': 1, 'header': True, 'seed': 0 }
pDecideK.lang    = "Rscript"
pDecideK.script  = """
set.seed ({{proc.args.seed}})
library ('factoextra')
data = read.table ("{{infile}}", row.names={{proc.args.rownames}}, header={{proc.args.header | Rbool }}, check.names=F)
data = data[, apply(data, 2, sd)!=0, drop=F]
k = 0
parseK = function (k) {
	parts = noquote (unlist(strsplit(k, ":")))
	return (c(as.numeric(parts[2]), as.numeric(parts[3])))
}
argsk = "{{proc.args.method}}"
if (startsWith (argsk, "elbow")) {
	library ('GMD')
	if (argsk == "elbow") argsk = "elbow:.95:.01"
	parts = parseK (argsk)
	d     = dist (data)
	cs    = css.hclust (d)
	elb   = elbow.batch(cs, ev.thres = parts[1], inc.thres = parts[2])
	k     = elb$k
} else if (startsWith (argsk, "pamk")) {
	library(fpc)
	if (argsk == 'pamk') argsk = "pamk:2:15"
	parts = parseK (argsk)
	p     = pamk (data, krange=parts[1]:parts[2])
	k     = p$nc
} else if (startsWith (argsk, "calinski")) {
	require(vegan)
	if (argsk == 'calinski') argsk = "calinski:2:15"
	parts = parseK (argsk)
	fit   = cascadeKM(scale(data, center = TRUE,  scale = TRUE), parts[1], parts[2], iter = 1000)
	k     = as.numeric(which.max(fit$results[2,]))
} else if (startsWith (argsk, "mclust")) {
	library(mclust)
	if (argsk == 'mclust') argsk = "mclust:2:15"
	parts = parseK (argsk)
	d_cl  = Mclust(as.matrix(data), G=parts[1]:parts[2])
	k     = dim(d_cl$z)[2]
} else if (startsWith (argsk, "ap")) {
	library(apcluster)
	apcl  = apcluster(negDistMat(r=2), data)
	k     = length(apcl@clusters)
} else if (startsWith (argsk, "gap")) {
	library (cluster)
	if (argsk == 'gap') argsk = "gap:2:15"
	parts = parseK (argsk)
	tabs  = clusGap(data, kmeans, parts[2], B = mean(10*parts[2], 500))
	tabs  = tabs$Tab
	for (i in parts[1]:(parts[2]-1)) {
		k = i
		if (tabs[i+1, "gap"] - tabs[i, "gap"] < 0) break
	}
} else if (startsWith (argsk, "nbclust")) {
	library(NbClust)
	if (argsk == 'nbclust') argsk = "nbclust:2:15"
	parts = parseK (argsk)
	nb    = NbClust(data, diss=NULL, distance = "euclidean", 
					min.nc=parts[1], max.nc=parts[2], method = "kmeans")
	nb    = nb$Best.nc[1, ]
	nb    = as.data.frame (table(nb))
	k     = nb[which.max (nb$Freq), 1]
} else {
	k     = as.numeric(argsk)
}
k         = min (nrow(data), k)
write (k, "{{kfile}}")
"""



"""
@name:
	pKMeans
@description:
	Do k-means clustering
@input:
	`infile:file`:    The input coordinates of the points.
	`k`:              Number of clusters, it could also be a file with the number in it.
@output:
	`outdir:dir`: The output of final results
@args:
	`rownames`:       The `row.names` for `read.table` to read the input file, default: 1.
	`header`:         The `header` argument for `read.table` to read the input file, default: True.
	`algorithm`:      The `algorithm` argument for `kmeans`, default "Hartigan-Wong" (could also be "Lloyd", "Forgy", "MacQueen")
	`niter`:          The `max.iter` argument for `kmeans`, default: 10.
	`nstart`:         The `nstart` argument for `kmeans`, default: 25.
	`caption`:        The caption for the `fviz_cluster`, default: "K-means with K=%k%".
@requires:
	[`r-factoextra`](https://cran.r-project.org/web/packages/factoextra/index.html)
"""
pKMeans = proc ()
pKMeans.input   = 'infile:file, k'
pKMeans.output  = 'outdir:dir:{{infile | fn}}.kmeans'
pKMeans.args    = { 'rownames': 1, 'header': True, 'algorithm': 'Hartigan-Wong', 'niter': 10, 'nstart': 25, 'caption': 'K-means with K=%K%' }
pKMeans.lang    = 'Rscript'
pKMeans.script  = """
library ('factoextra')
data      = read.table ("{{infile}}", row.names={{proc.args.rownames}}, header={{proc.args.header | Rbool }}, check.names=F)
data      = data[, apply(data, 2, sd)!=0, drop=F]
k         = as.numeric ("{{k}}")
if (is.na (k)) {
	kf    = "{{k}}"
	k     = as.numeric (readChar (kf, file.info(kf)$size))
}
km        = kmeans (data, k, iter.max={{proc.args.niter}}, nstart="{{proc.args.nstart}}", algorithm="{{proc.args.algorithm}}")
clustfile = file.path ("{{outdir}}", "cluster.txt")
write.table (km$cluster, clustfile, quote=F, col.names = F, sep="\\t")
clustfig  = file.path ("{{outdir}}", "cluster.png")
png (file = clustfig)
print(fviz_cluster (km, data=data, main = gsub("%K%", k, "{{proc.args.caption}}")))
dev.off()
"""

"""
@name:
	pPamk
@description:
	Do clustering using [fpc::pamk](https://www.rdocumentation.org/packages/fpc/versions/2.1-10/topics/pamk)
@input:
	`infile:file`:  The input coordinate file
@output:
	`outdir:dir`:   The output directory
@args:
	`rownames`:     The `row.names` for `read.table` to read the input file, default: 1.
	`header`:       The `header` argument for `read.table` to read the input file, default: True.
	`min`:          The min # clusters to try, default: 2
	`max`:          The max # clusters to try, default: 15
	`caption`:      The caption for the `fviz_cluster`, default: "Partitioning Around Medoids (K=%K%)".
	`seed`:         The seed for randomization, default: 0.
@requires:
	[`r-factoextra`](https://cran.r-project.org/web/packages/factoextra/index.html)
	[`r-fpc`](https://cran.r-project.org/web/packages/fpc/index.html)
"""
pPamk = proc()
pPamk.input    = "infile:file"
pPamk.output   = "outdir:dir:{{infile | fn}}.pamk"
pPamk.args     = {'rownames': 1, 'header': True, "min":2, "max":15, "seed":0, "caption": "Partitioning Around Medoids (K=%k%)"}
pPamk.lang     = "Rscript"
pPamk.script   = """
library(fpc)
library(factoextra)
data = read.table ("{{infile}}", row.names={{proc.args.rownames}}, header={{proc.args.header | Rbool }}, check.names=F)
data = data[, apply(data, 2, sd)!=0, drop=F]
p = pamk (data, krange={{proc.args.min}}:{{proc.args.max}}, seed={{proc.args.seed}})
p$pamobject$data = data
clustfile = file.path ("{{outdir}}", "cluster.txt")
write.table (p$pamobject$clustering, clustfile, quote=F, col.names = F, sep="\\t")
clustfig  = file.path ("{{outdir}}", "cluster.png")
png (file = clustfig)
print (fviz_cluster (p$pamobject, main = gsub("%K%", p$nc, "{{proc.args.caption}}")))
dev.off()
"""

"""
@name:
	pClara
@description:
	CLARA is a partitioning method used to deal with much larger data sets (more than several thousand observations) in order to reduce computing time and RAM storage problem.
@input:
	`infile:file`:  The input coordinate file
	`k`:            Number of clusters, it could also be a file with the number in it.
@output:
	`outdir:dir`:   The output of final results
@args:
	`rownames`:     The `row.names` for `read.table` to read the input file, default: 1.
	`header`:       The `header` argument for `read.table` to read the input file, default: True.
	`samples`:      The `samples` argument for `clara`, default: 5.
	`caption`:      The caption for the `fviz_cluster`, default: "CLARA Clustering with K=%k%".
@requires:
	[`r-cluster`](https://cran.r-project.org/web/packages/cluster/index.html)
	[`r-factoextra`](https://cran.r-project.org/web/packages/factoextra/index.html)
"""
pClara  = proc()
pClara.input   = "infile:file, k"
pClara.output  = 'outdir:dir:{{infile | fn}}.clara'
pClara.args    = { 'rownames': 1, 'header': True, 'samples': 5, 'caption': 'CLARA Clustering with K=%K%' }
pClara.lang    = 'Rscript'
pClara.script  = """
library ('factoextra')
library ('cluster')
data = read.table ("{{infile}}", row.names={{proc.args.rownames}}, header={{proc.args.header | Rbool }}, check.names=F)
data = data[, apply(data, 2, sd)!=0, drop=F]
k         = as.numeric ("{{k}}")
if (is.na (k)) {
	kf    = "{{k}}"
	k     = as.numeric (readChar (kf, file.info(kf)$size))
}

clr       = clara (data, k=k, samples={{proc.args.samples}})
clustfile = file.path ("{{outdir}}", "cluster.txt")
write.table (clr$cluster, clustfile, quote=F, col.names = F, sep="\\t")
clustfig  = file.path ("{{outdir}}", "cluster.png")
png (file = clustfig)
labels    = colnames(data)
print (fviz_cluster (clr, main = gsub("%K%", k, "{{proc.args.caption}}")))
dev.off()
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
data = read.table ("{{infile}}", row.names={{proc.args.rownames}}, header={{proc.args.header | Rbool }}, check.names=F)
data = data[, apply(data, 2, sd)!=0, drop=F]
mc        = Mclust (data, G={{proc.args.min}}:{{proc.args.max}})
clustfile = file.path ("{{outdir}}", "cluster.txt")
write.table (mc$classification, clustfile, quote=F, col.names = F, sep="\\t")
clustfig  = file.path ("{{outdir}}", "cluster.png")
png (file = clustfig)
labels    = colnames(data)
print (fviz_cluster (mc, main = gsub("%K%", max(mc$classification), "{{proc.args.caption}}")))
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
data = read.table ("{{infile}}", row.names={{proc.args.rownames}}, header={{proc.args.header | Rbool }}, check.names=F)
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
print (fviz_cluster (fvizobj, main = gsub("%K%", nclust, "{{proc.args.caption}}")))
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
	[`r-fastcluster`](https://cran.r-project.org/web/packages/fastcluster/index.html) if `proc.args.fast` is True
	[`r-ggdendro`](https://cran.r-project.org/web/packages/ggdendro/index.html) if `proc.args.gg` is True
"""
pHClust = proc()
pHClust.input  = "infile:file"
pHClust.output = "outdir:dir:{{infile | fn}}.hclust"
pHClust.args   = {"fast":False, "gg":False, "rownames":1, "header":True, 'method': 'complete', 'rotate': False, 'transpose': False}
pHClust.lang   = "Rscript"
pHClust.script = """
data = read.table ("{{infile}}", row.names={{proc.args.rownames}}, header={{proc.args.header | Rbool }}, check.names=F)
if ({{proc.args.transpose | Rbool}}) data = t(data)
data = data[, apply(data, 2, sd)!=0, drop=F]
dmat = dist(data)
if ({{proc.args.fast | Rbool}}) {
	library('fastcluster')
}
orderfile = file.path ("{{outdir}}", "hclust.order.txt")
mergefile = file.path ("{{outdir}}", "hclust.merge.txt")
clustfig  = file.path ("{{outdir}}", "hclust.png")
hobj = hclust (dmat, method="{{proc.args.method}}")
png (file = clustfig)
if ({{proc.args.gg | Rbool}}) {
	library('ggplot2')
	library('ggdendro')
	ggdendrogram (hobj, rotate={{proc.args.rotate | Rbool}})
} else {
	plot (as.dendrogram(hobj), horiz={{proc.args.rotate | Rbool}})
}
dev.off()
orderobj = data.frame (Order=hobj$order, Labels=hobj$labels)
mergeobj = data.frame (Merge1=hobj$merge[,1], Merge2=hobj$merge[,2], Height=hobj$height)
write.table (orderobj, orderfile, quote=F, row.names=F, col.names = T, sep="\\t")
write.table (mergeobj, mergefile, quote=F, row.names=F, col.names = T, sep="\\t")

"""
