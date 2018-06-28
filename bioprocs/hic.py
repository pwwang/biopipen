"""
Some long range interaction processes, not necessarily Hi-C
"""
from pyppl import Proc
from . import params

"""
@name:
	pPartners
@description:
	Find the interaction partners of the regions in input file.
@input:
	`regfile:file`: The region file for the regions to find partners for.
	`intfile:file`: The interaction file
@output:
	`outfile:file`: The regions with partners.
@args:
	`regtype`: The type of region file. Default: `auto` (tell from file extension)
		- Could also be `bed` or `bedx`
	`inttype`: The type of interaction file. Default: `auto`
		- Could also be `bedpe`, `chiapet.tool`, `hiclib` and `bed12`
"""
pPartners                = Proc(desc = 'Find the interaction partners of the regions in input file.')
pPartners.input          = "regfile:file, intfile:file"
pPartners.output         = "outfile:file:{{in.regfile | fn2}}-{{in.intfile | fn2}}.partners.bedx"
pPartners.args.regtype   = "auto" # bed,   bedx
pPartners.args.inttype   = "auto" # bedpe, chiapet.tool, hiclib, bed12
pPartners.lang           = params.python.value
pPartners.script         = "file:scripts/hic/pPartners.py"