from pyppl import Box
from bioprocs.utils.tsvio import tsvops

infile  = {{in.infile | quote}}
outfile = {{out.outfile | quote}}
inopts  = {{args.inopts}}
outopts = {{args.outopts}}

opshelper = [line for line in {{args.opshelper | quote}}.splitlines() if line]
while opshelper and (all([line and line[0] == ' ' for line in opshelper]) or all([line and line[0] == '\t' for line in opshelper])):
	opshelper = [line[1:] for line in opshelper]
exec('\n'.join(opshelper))

ops = {{args.ops}}

tsvops(infile, outfile, inopts, outopts, ops)

