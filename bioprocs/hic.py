"""
Some long range interaction processes, not necessarily Hi-C
"""
from pyppl import Proc, Box
from . import params
from . import delefactory, procfactory
from modkit import Modkit
Modkit().delegate(delefactory())

@procfactory
def _pPartners():
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
	pPartners.output         = "outfile:file:{{i.regfile | fn2}}-{{i.intfile | fn2}}.partners.bedx"
	pPartners.args.regopts   = Box(ftype = "auto") # bed,   bedx
	pPartners.args.intopts   = Box(ftype = "auto") # bedpe, chiapet.tool, hiclib, bed12
	pPartners.lang           = params.python.value
	pPartners.script         = "file:scripts/hic/pPartners.py"
	return pPartners

