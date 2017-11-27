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
pRWR.output     = "outfile:file:{{in.Wfile | fn}}-{{in.Efile | fn}}.rwr"
pRWR.args.c     = 0.1
pRWR.args.eps   = 1e-5
pRWR.args.niter = 10000
pRWR.args.normW = True
pRWR.args.normE = True
pRWR.lang       = params.Rscript.value
pRWR.script     = "file:scripts/algorithm/pRWR.r"




