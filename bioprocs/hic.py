"""
Some long range interaction processes, not necessarily Hi-C
"""
from pyppl import Proc
from . import params
from .utils import read, write

pPartners                = Proc(desc = 'Find the interaction partners of the regions in input file.')
pPartners.input          = "regfile:file, intfile:file"
pPartners.output         = "outfile:file:{{in.regfile | fn2}}-{{in.intfile | fn2}}.partners.bedx"
pPartners.args.regtype   = "auto" # bed,   bedx
pPartners.args.inttype   = "auto" # bedpe, chiapet.tool, hiclib, bed12
pPartners.envs.readBed12 = read.bed12.py
pPartners.envs.readBed   = read.bed.py
pPartners.envs.readBedpe = read.bedpe.py
pPartners.envs.readBedx  = read.bedx.py
pPartners.envs.writeBedx = write.bedx.py
pPartners.lang           = params.python.value
pPartners.script         = "file:scripts/hic/pPartners.py"