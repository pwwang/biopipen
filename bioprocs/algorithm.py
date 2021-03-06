"""A set of algorithms or models"""
from diot import Diot
from pyppl import Proc
from . import params, proc_factory

pRWR = proc_factory(
	desc   = 'Do random walk with restart (RWR).',
	config = Diot(annotate = """
		@input:
			`Wfile:file`: The adjecent matrix
			`Efile:file`: The start vector
		@output:
			`outfile:file`: The output of final probabilities
		@args:
			`c`:       The restart probability. Default: 0.1
			`eps`:     The convergent cutoff || R(i+1) - R(i) ||. Default: 1e-5
			`niter`:   Max iterations to stop. Default: 10000
			`normW`:   Weather to normalize W or not, default True.
				- Laplacian normalization is used (more to add).
			`normE`:   Weather to normalize E or not, default True.
				- E will be normalized as: E = E/sum(E)
		@requires:
			[NetPreProc]
			desc = "Package for the pre-processing and normalization of graphs."
			url = "https://cran.r-project.org/web/packages/NetPreProc/index.html"
			when = "[[ {{args.normW}} == True ]]"
			validate = '[[ $({{proc.lang}} --vanilla -e 'library(NetPreProc)' 2>&1) == \
				*"Network Pre-Processing package"* ]]'
			install = '{{proc.lang}} -e \'install.packages("NetPreProc", \
				repos="https://cran.rstudio.com")\''
		"""))
pRWR.input      = "Wfile:file, Efile:file"
pRWR.output     = "outfile:file:{{i.Wfile | fn2}}.rwr.txt"
pRWR.lang       = params.Rscript.value
pRWR.args.c     = 0.1
pRWR.args.eps   = 1e-5
pRWR.args.niter = 10000
pRWR.args.normW = True
pRWR.args.normE = True

pAR = Proc(
	desc   = 'Affinity Regression.',
	config = Diot(annotate = """
	@name:
		pAR
	@description:
		Affinity Regression.
		Ref: https://www.nature.com/articles/nbt.3343
		```
				b           c        d          d
			_________    _______    ____       ____
			|       |    |  W  |    |  |       |  |
			a |   D   |  b |_____|  c |Pt|  =  a |Y |   <=>
			|_______|               |__|       |  |
											|__|
		kronecker(P, YtD)*vec(W) = vec(YtY)             <=>
		X*vec(W) = vec(YtY)
		WPt:
			c        d              d
			_______    ____          _____
			|  W  |    |  |          |   |
			b |_____|  c |Pt|  --->  b |___|
					|__|
		YtDW:
		WtDtY:
			b           a        d               d
			_______    _________   ____           _____
			|  Wt |    |       |   |  |           |   |
		c |_____|  b |   Dt  | a |Y |    ---> c |___|
					|_______|   |  |
								|__|
		```
	@input:
		`D:file` : The D matrix
		`Pt:file`: The Pt matrix
		`Y:file`:  The Y matrix
			- All input files could be gzipped
	@output:
		`W:file`:  The interaction matrix
		`outdir:dir`: The output directory
	@args:
		`seed`:  The seed for sampling the training set.
		`tfrac`: The fraction of samples used for training.
	"""))
pAR.input  = 'D:file, Pt:file, Y:file'
pAR.output = [
	'W:file:{{i.D | fn}}-{{i.Pt | fn}}-{{i.Y | fn}}.AR/W.txt',
	'outdir:dir:{{i.D | fn}}-{{i.Pt | fn}}-{{i.Y | fn}}.AR'
]
pAR.lang         = params.Rscript.value
pAR.args.seed    = None
pAR.args.tfrac   = .5
pAR.args.inopts  = Diot(cnames = True, rnames = True)
pAR.args.svdP    = 0
pAR.args.predY   = True
pAR.args.WPt     = True
pAR.args.WtDtY   = True
pAR.args.nfold   = 3
pAR.args.nthread = 1
pAR.args.method  = 'glmnet' # admm

pColoc = Proc(
	desc   = "Bayes Factor colocalisation analyses using R `coloc` package.",
	config = Diot(annotate = """
	@description:
		Bayes Factor colocalisation analyses using R `coloc` package.
		`coloc` package can accept multiple formats of input. Here we adopt the one using pvalues.
		`coloc.abf(dataset1=list(pvalues=p1,N=nrow(X1),type="quant"), dataset2=list(pvalues=p2,N=nrow(X2),type="quant"), MAF=maf)`
	@input:
		`infile:file`: The input file including the MAF, pvalues of 1st and 2nd phenotypes
			- The first 6 columns are in BED6 format.
			- 7th : MAF
			- 8th : Pvalues for the 1st phenotype
			- 9th : Pvalues for the 2nd phenotype
			- This file could have a header with the names for phenotypes
			- Snps have to be on the same chromosome, and sorted by positions.
	@output:
		`outfile:file`: The output file including:
			- # snps, PP.H0.abf, PP.H1.abf, PP.H2.abf, PP.H3.abf and PP.H4.abf
		`outdir:dir`  : The output directory containing the output file and plots.
	@args:
		`plot`: Do manhattan plot? Default: `True`
	"""))
pColoc.input  = 'infile:file'
pColoc.output = [
	'outfile:file:{{i.infile | fn2}}.coloc/{{i.infile | fn2}}.coloc.txt',
	'outdir:dir:{{i.infile | fn2}}.coloc'
]
pColoc.args.inopts  = Diot(cnames = True, rnames = False)
pColoc.args.plot    = True
pColoc.args.ggs     = Diot()
pColoc.args.params  = Diot()
pColoc.args.devpars = Diot(res = 300, height = 2000, width = 2000)
pColoc.args.hifile  = ''
pColoc.args.hilabel = True
pColoc.lang         = params.Rscript.value
