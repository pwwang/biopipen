from os import path
from pyppl import Box
from gff import Gff
from bioprocs.utils.tsvio2 import TsvWriter, TsvReader, TsvRecord
from bioprocs.utils import logger

infile   = {{ i.infile | quote}}
outfile  = {{ o.outfile | quote}}
notfound = {{ args.notfound | quote}}
genecol  = {{ args.genecol or 0 | repr}}
inopts   = {{ args.inopts | repr}}
refgene  = {{ args.refgene | quote}}

if not path.isfile(refgene):
	raise OSError('Refgene file does not exists: {}'.format(refgene))

# get genes
genes  = TsvReader(infile, **inopts).dump(genecol)
genes  = dict(zip(genes, [False] * len(genes)))
writer = TsvWriter(outfile)
writer.cnames = ['CHR', 'START', 'END', 'NAME', 'SCORE', 'STRAND']

gff    = Gff(refgene)
for g in gff:
	attrs = g['attributes']
	if attrs['gene_id'] not in genes:
		continue
	r        = TsvRecord()
	r.CHR    = g['seqid']
	r.START  = g['start']
	r.END    = g['end']
	r.SCORE  = g['score']
	r.STRAND = g['strand']
	r.NAME   = attrs['gene_id']
	writer.write(r)
writer.close()

for g, v in genes.items():
	if not v:
		logger.warning('Gene: {!r} not found.'.format(g))
