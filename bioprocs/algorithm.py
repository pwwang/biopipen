from pyppl import Proc
from . import params
"""
@name:
	pRWR
@description:
	Do random walk with restart (RWR)
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
	if normW = True, R package `NetPreProc` is required.
"""
pRWR            = Proc (desc = 'Do random walk with restart (RWR).')
pRWR.input      = "Wfile:file, Efile:file"
pRWR.output     = "outfile:file:{{in.Wfile | fn2}}.rwr.txt"
pRWR.args.c     = 0.1
pRWR.args.eps   = 1e-5
pRWR.args.niter = 10000
pRWR.args.normW = True
pRWR.args.normE = True
pRWR.lang       = params.Rscript.value
pRWR.script     = "file:scripts/algorithm/pRWR.r"


"""
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
	``
"""
pAR            = Proc(desc =  'Affinity Regression.')
pAR.input      = 'D:file, Pt:file, Y:file'
pAR.output     = [
	'W:file:{{in.D | fn}}-{{in.Pt | fn}}-{{in.Y | fn}}.AR/W.txt', 
	'outdir:dir:{{in.D | fn}}-{{in.Pt | fn}}-{{in.Y | fn}}.AR'
]
pAR.args.seed     = -1
pAR.args.tfrac    = .5
pAR.args.normD    = 'L2'
pAR.args.normPt   = 'meanCenter'
pAR.args.normY    = 'meanCenter'
pAR.args.svdP     = 0
pAR.args.predY    = True
pAR.args.WPt      = True
pAR.args.WtDtY    = True
pAR.args.nfolds   = 3
pAR.args.parallel = False
pAR.args.method   = 'glmnet' # admm
pAR.lang          = params.Rscript.value
pAR.script        = "file:scripts/algorithm/pAR.r"





