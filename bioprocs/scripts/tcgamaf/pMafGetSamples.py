from pyppl import Box
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter

infile  = {{ i.infile | quote}}
outfile = {{o.outfile | quote}}
out     = {{args.out | quote}}
outopts = {{args.outopts |repr}}

CNAMES = Box(
	tumor  = ['TUMOR'],
	normal = ['NORMAL'],
	both   = ['TUMOR', 'NORMAL']
)


reader = TsvReader(infile)
writer = TsvWriter(outfile)
cnames = outopts.get('cnames')
if cnames is True:
	cnames = CNAMES[out]
if cnames:
	writer.cnames = cnames
	writer.writeHead()

prev_tumor = prev_normal = None
for r in reader:
	tumor, normal = r.Tumor_Sample_Barcode, r.Matched_Norm_Sample_Barcode
	if prev_tumor == tumor and prev_normal == normal:
		continue
	if out == 'tumor':
		writer.write([tumor])
	elif out == 'normal':
		writer.write([normal])
	else:
		writer.write([tumor, normal])
	prev_tumor = tumor
	prev_normal = normal
