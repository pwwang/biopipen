from pyppl import Proc, Box
from bioprocs import params, rimport

"""
@name:
	pGEP70
@description:
	Calculate GEP70 scores for multiple mylenoma 70-gene-signatures
"""
pGEP70              = Proc(desc = 'Calculate GEP70 scores for multiple mylenoma 70-gene-signatures')
pGEP70.input        = 'infile:file' # make sure the expression values are log2-scale normalized
pGEP70.output       = 'outfile:file:{{in.infile | fn2}}.gep70.txt'
pGEP70.args.gep70   = params.gep70.value
pGEP70.args.inopts  = Box(cnames = True)
pGEP70.args.lang    = 'Rscript'
pGEP70.envs.rimport = rimport
pGEP70.script       = "file:scripts/misc/pGEP70.r"