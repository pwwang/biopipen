"""Processes for clustering data"""
from pyppl import Proc, Box
from . import params, rimport
from . import delefactory, procfactory
from modkit import Modkit
Modkit().delegate(delefactory())

@procfactory
def _pDist2Feats():
	"""
	@name:
		pDist2Feats
	@description:
		Convert a distance matrix to coordinates, using multidimensional scaling.
	@input:
		`infile:file`: The distance matrix, could be a full distance matrix, a triangle matrix or a pair-wise distance file
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
	@output:
		`outfile:file`: The output coordinate file
	@args:
		`informat`: The format of the input file: full, triangle or pair. Default: full
			- Could also be upper, lower, pair
		`k`:        How many dimensions? Default: 2 (R^2)
	"""
	pDist2Feats            = Proc(desc = 'Convert a distance matrix to coordinates, using multidimensional scaling.')
	pDist2Feats.input      = "infile:file"
	pDist2Feats.output     = "outfile:file:{{i.infile | fn}}.xyz"
	pDist2Feats.args.infmt = 'full'
	pDist2Feats.args.k     = 2
	pDist2Feats.lang       = params.Rscript.value
	pDist2Feats.script     = "file:scripts/cluster/pDist2Feats.r"
	return pDist2Feats

@procfactory
def _pFeats2Dist():
	"""
	@name:
		pFeats2Dist
	@description:
		Calculate the distance of pairs of instances
	@input:
		`infile:file`: The input file. Rows are instances, columns are features.
	@output:
		`outfile:file`: The output distance matrix
	@args:
		`inopts`: the options for input file
			- `cnames`: Whether input file has column names
			- `rnames`: Whether input file has row names
		`transfm`: Transform the results in follow way. Default: `scale`
			- `scale`: scale the distance to [0, 1]
			- `similarity`: 1 - scaled distance
			- `False`: don't do any transformation.
		`na` : Replace NAs with? Default: `0`
		`method`: Method used to calculate the distance. See available below. Default: `euclidean`
			- euclidean, manhattan, minkowski, chebyshev, sorensen, gower, soergel, kulczynski_d,
			- canberra, lorentzian, intersection, non-intersection, wavehedges, czekanowski, motyka,
			- kulczynski_s, tanimoto, ruzicka, inner_product, harmonic_mean, cosine, hassebrook,
			- jaccard, dice, fidelity, bhattacharyya, hellinger, matusita, squared_chord,
			- quared_euclidean, pearson, neyman, squared_chi, prob_symm, divergence, clark,
			- additive_symm, kullback-leibler, jeffreys, k_divergence, topsoe, jensen-shannon,
			- jensen_difference, taneja, kumar-johnson, avg
			- see: https://cran.r-project.org/web/packages/philentropy/vignettes/Distances.html
	@requires:
		`r-philentropy`
	"""
	pFeats2Dist              = Proc(desc = 'Calculate the distance between each pair of rows')
	pFeats2Dist.input        = 'infile:file'
	pFeats2Dist.output       = 'outfile:file:{{i.infile | fn}}.dist.txt'
	pFeats2Dist.args.inopts  = Box(cnames = True, rnames = True)
	pFeats2Dist.args.transfm = 'scale'
	pFeats2Dist.args.na      = 0
	pFeats2Dist.args.method  = 'euclidean'
	pFeats2Dist.envs.rimport = rimport
	pFeats2Dist.lang         = params.Rscript.str()
	pFeats2Dist.script       = "file:scripts/cluster/pFeats2Dist.r"
	return pFeats2Dist

@procfactory
def _pCluster():
	"""
	@name:
		pCluster
	@description:
		Use `optCluster` to select the best number of clusters and cluster method, then perform the clustering
	@input:
		`infile:file`: The input matrix file. Clustering will be performed against rows. If not, set `args.transpose` = True
	@output:
		`outfile:file`: The output cluster file
		`outdir:dir`:   The output directory containing the figures
	@args:
		`transpose`:    Transpose the input matrix. Default: False
		`cnames`:       Whether the input matrix contains header before transposing. Default: False
		`rnames`:       Which column is the rownames before transposing. Default: 1
		`plot`:         Whether plot the cluster. Default: True
		`minc`:         Min number of clusters to test. Default: 2
		`maxc`:         Min number of clusters to test. Default: 15
			- If number of rows (nrows) <= 15, then max = nrows - 1
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
	pCluster                = Proc (desc = 'Use `optCluster` to select the best number of clusters and cluster method, then perform the clustering.')
	pCluster.input          = 'infile:file'
	pCluster.output         = 'outfile:file:{{i.infile | fn}}.cluster/clusters.txt, outdir:dir:{{i.infile | fn}}.cluster'
	pCluster.args.transpose = False
	pCluster.args.cnames    = True
	pCluster.args.rnames    = True
	pCluster.args.plot      = True
	pCluster.args.minc      = 2
	pCluster.args.maxc      = 15
	# "agnes", "clara", "diana", "fanny", "hierarchical", "kmeans", "model", "pam", "som", "sota", "em.nbinom", "da.nbinom", "sa.nbinom", "em.poisson", "da.poisson", "sa.poisson"
	# Note fanny reports error and the algorithm cannot go further.
	# It is excluded by default, you can add it manually
	# model and sota are excluded because they are slow
	# You can also manually add them
	pCluster.args.methods   = 'all'
	pCluster.args.isCount   = False
	pCluster.args.devpars   = Box(res=300, width=2000, height=2000)
	pCluster.lang           = params.Rscript.value
	pCluster.script         = "file:scripts/cluster/pCluster.r"
	return pCluster


@procfactory
def _pMCluster():
	"""
	@name:
		pMCluster
	@description:
		Use `r-mclust` to do clustering. Current just do simple clustering with the package
	@input:
		`infile:file`: The input a coordinate file
	@output:
		`outdir:dir`:  The output of final results
	@args:
		`transpose`: Transpose the input matrix? Default: False
		`rnames`:  The `row.names` for `read.table` to read the input file, default: True.
		`cnames`:  The `header` argument for `read.table` to read the input file, default: True.
		`caption`: The caption for the `fviz_cluster`, default: "CLARA Clustering".
		`minc`:    The min # clusters to try, default: 2
		`maxc`:    The max # clusters to try, default: 15
	@requires:
		[`r-mclust`](https://cran.r-project.org/web/packages/mclust/index.html)
		[`r-factoextra`](https://cran.r-project.org/web/packages/factoextra/index.html)
	"""
	pMCluster                = Proc(desc = 'Use `r-mclust` to do clustering. Current just do simple clustering with the package.')
	pMCluster.input          = "infile:file"
	pMCluster.output         = 'outfile:file:{{i.infile | fn}}.mclust/clusters.txt, outdir:dir:{{i.infile | fn}}.mclust'
	pMCluster.args.transpose = False
	pMCluster.args.rnames    = True
	pMCluster.args.cnames    = True
	pMCluster.args.minc      = 2
	pMCluster.args.maxc      = 15
	pMCluster.args.devpars   = Box(res=300, width=2000, height=2000)
	pMCluster.lang           = params.Rscript.value
	pMCluster.script         = "file:scripts/cluster/pMCluster.r"
	return pMCluster


@procfactory
def _pAPCluster():
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
		`transpose`: Transpose the input matrix? Default: False
		`rownames`:     The `row.names` for `read.table` to read the input file, default: 1.
		`header`:       The `header` argument for `read.table` to read the input file, default: True.
		`caption`:      The caption for the `fviz_cluster`, default: "APClustering".
	@requires:
		[`r-apcluster`](https://cran.r-project.org/web/packages/apcluster/index.html)
		[`r-factoextra`](https://cran.r-project.org/web/packages/factoextra/index.html)
	"""
	pAPCluster                = Proc(desc = 'Use `r-apcluster` to do clustering.')
	pAPCluster.input          = "infile:file"
	pAPCluster.output         = 'outfile:file:{{i.infile | fn}}.apclust/clusters.txt, outdir:dir:{{i.infile | fn}}.apclust'
	pAPCluster.args.transpose = False
	pAPCluster.args.rnames    = True
	pAPCluster.args.cnames    = True
	pAPCluster.args.devpars   = Box(res=300, width=2000, height=2000)
	pAPCluster.lang           = params.Rscript.value
	pAPCluster.script         = "file:scripts/cluster/pAPCluster.r"
	return pAPCluster

@procfactory
def _pHCluster():
	"""
	@name:
		pHCluster
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
	pHCluster = Proc(desc = 'Do hierarchical clustering.')
	pHCluster.input  = "infile:file"
	pHCluster.output = [
		"outmerge:file:{{i.infile | fn}}.hclust/merge.txt",
		"outorder:file:{{i.infile | fn}}.hclust/order.txt",
		"outdir:dir:{{i.infile | fn}}.hclust"
	]
	pHCluster.args.transpose = False
	pHCluster.args.cnames    = True
	pHCluster.args.rnames    = True
	pHCluster.args.method    = 'complete'
	pHCluster.args.rotate    = False
	pHCluster.args.fast      = False
	pHCluster.args.devpars   = Box(res=300, width=2000, height=2000)
	pHCluster.lang           = params.Rscript.value
	pHCluster.script         = "file:scripts/cluster/pHCluster.r"
	return pHCluster

@procfactory
def _pKMeans():
	"""
	@name:
		pKMeans
	"""
	pKMeans        = Proc(desc = 'Do K-Means clustering')
	pKMeans.input  = 'infile:file'
	pKMeans.output = [
		'outfile:file:{{i.infile | fn}}.kmeans/{{i.infile | fn}}.kmeans.txt',
		'outdir:dir:{{i.infile | fn}}.kmeans'
	]
	pKMeans.args.inopts  = Box(cnames = True, rnames = True)
	pKMeans.args.k       = None
	pKMeans.args.plot    = True
	pKMeans.args.ggs     = Box()
	pKMeans.args.devpars = Box(res=300, width=2000, height=2000)
	pKMeans.envs.rimport = rimport
	pKMeans.lang         = params.Rscript.str()
	pKMeans.script       = "file:scripts/cluster/pKMeans.r"
	return pKMeans

