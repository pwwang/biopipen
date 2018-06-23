from pyppl import Proc, Box
from bioprocs import params, rimport

"""
@name:
	pGEP70
@description:
	Calculate GEP70 scores for multiple mylenoma 70-gene-signatures
@input:
	`infile:file`: The input file with expression matrix
		- Columns are samples, rows are genes
@output:
	`outfile:file`: The output files with gep70 scores for each sample.
		- Samples become rows, just one column is in the file.
@args:
	`inopts`: The input options.
		- `cnames`: Whether the input file has column names. Default: `True`
	`gep70`: The GEP70 genes. 
		- Column 1: up-regulated genes (51)
		- Column 2: down-regulated genes (19)
"""
pGEP70              = Proc(desc = 'Calculate GEP70 scores for multiple mylenoma 70-gene-signatures')
pGEP70.input        = 'infile:file' # make sure the expression values are log2-scale normalized
pGEP70.output       = 'outfile:file:{{in.infile | fn2}}.gep70.txt'
pGEP70.args.gep70   = params.gep70.value
pGEP70.args.inopts  = Box(cnames = True)
pGEP70.args.lang    = 'Rscript'
pGEP70.envs.rimport = rimport
pGEP70.script       = "file:scripts/misc/pGEP70.r"