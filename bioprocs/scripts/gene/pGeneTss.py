import sys

from pyppl import Box
from bioprocs.utils.tsvio import TsvWriter, TsvRecord
from bioprocs.utils.gene import genenorm
from bioprocs.utils import logger

genome   = {{args.genome | quote}}
poscol   = 'genomic_pos_%s' % genome if genome!= 'hg38' else 'genomic_pos'
notfound = {{args.notfound | quote}}
# get the genes
genes = genenorm(
	{{i.infile | quote}},
	notfound = notfound,
	frm      = {{args.frm | quote}},
	to       = "%s,symbol" % poscol,
	genome   = genome,
	cachedir = {{args.cachedir | quote}},
	inopts   = {{args.inopts}},
	genecol  = {{args.genecol | quote}}
)
outopts = {{args.outopts}}
writer = TsvWriter({{o.outfile | quote}}, **outopts)
if outopts['query']:
	writer.meta.add('QUERY')
if outopts['head']:
	writer.writeHead(delimit = outopts['headDelimit'], prefix = outopts['headPrefix'], transform = outopts['headTransform'])
for gene, hit in genes.items():
	if not poscol in hit or not hit[poscol]:
		if notfound == 'skip':
			logger.warn('Gene not found: %s' % gene)
			continue
		else:
			raise ValueError('Gene not found: %s' % gene)

	pos      = hit[poscol]
	if isinstance(pos, list): pos = pos[0]
	r        = TsvRecord()
	r.CHR    = 'chr' + str(pos['chr'])
	r.START  = pos['start'] if pos['strand'] == 1 else pos['end']
	r.END    = (r.START + 1) if pos['strand'] == 1 else int(pos['end']) + 1
	r.NAME   = hit['symbol']
	r.STRAND = '+' if pos['strand'] == 1 else '-'
	r.SCORE  = '0'
	r.QUERY  = gene
	writer.write(r)
