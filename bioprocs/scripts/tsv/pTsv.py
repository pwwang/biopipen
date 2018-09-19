from pyppl import Box
from bioprocs.utils.tsvio import tsvops, TsvRecord

infile    = {{i.infile | quote}}
outfile   = {{o.outfile | quote}}
inopts    = {{args.inopts}}
outopts   = {{args.outopts}}
opshelper = {{args.opshelper | repr}}
if not isinstance(opshelper, list):
	opshelper = [opshelper]

opshelper = [line for line in opshelper if line]
exec('\n'.join(opshelper), globals())

ops = {{args.ops}}

tsvops(infile, outfile, inopts, outopts, ops)

