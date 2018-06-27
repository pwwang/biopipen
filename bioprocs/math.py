from pyppl import Proc, Box
from . import params, rimport

"""
@name:
	pRank
@description:
	Convert values to ranks.
@input:
	`infile:file`: The input file
@output:
	`outfile:file`: The output file with ranks.
@args:
	`na`: Where to put the `NA` values.
		- `"first"` : Put `NA` first
		- `"last"`  : Put `NA` last (default)
		- `"remove"`: Remove `NA` values
		- `"keep"`  : keep `NA` values
	`tie`: How to deal with ties
		- `"average"` : Use average ranks (default)
		- `"first"`   : Use the ranks come first
		- `"last"`    : Use the ranks come last
		- `"random"`  : Use the random ranks
		- `"max"`     : Use the max ranks
		- `"min"`     : Use the min ranks
	`byrow`: Calculate ranks by row (instead of by column)? Default: `True`
	`reverse`: Take the reverse rank? Default: `True`
		- Large number gets higher rank (smaller rank index)
		- `args.na` remains the same.
	`inopts`: The input options:
		- `cnames`: Whether the input file has header. Default: `True`
		- `rnames`: Whether the input file has row names. Default: `True`
		- `delimit`: The separator of columns. Default: `\t`
"""
pRank              = Proc(desc = 'Convert values to ranks')
pRank.input        = 'infile:file'
pRank.output       = 'outfile:file:{{in.infile | fn}}.rank.txt'
pRank.args.na      = 'last' # keep,         first,   remove
pRank.args.tie     = 'average' # "average", "first", "last", "random", "max", "min"
pRank.args.byrow   = True # else by column
pRank.args.reverse = True # large number ranks higher
pRank.args.inopts  = Box(
	cnames  = True,
	rnames  = True,
	delimit = "\t"
)
pRank.envs.rimport = rimport
pRank.lang         = params.Rscript.value
pRank.script       = "file:scripts/math/pRank.r"