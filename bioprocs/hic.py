"""
Some long range interaction processes, not necessarily Hi-C
"""
from pyppl import Proc
from . import params

pPartners                = Proc(desc = 'Find the interaction partners of the regions in input file.')
pPartners.input          = "regfile:file, intfile:file"
pPartners.output         = "outfile:file:{{in.regfile | fn2}}-{{in.intfile | fn2}}.partners.bedx"
pPartners.args.regtype   = "auto" # bed,   bedx
pPartners.args.inttype   = "auto" # bedpe, chiapet.tool, hiclib, bed12
pPartners.lang           = params.python.value
pPartners.script         = "file:scripts/hic/pPartners.py"