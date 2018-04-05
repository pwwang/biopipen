import sys, json

from bioprocs.utils.tsvio import TsvWriter, TsvReader, TsvRecord
from bioprocs.utils.gene import genenorm

infile  = {{in.infile   | quote}}
outfile = {{out.outfile | quote}}
inopts  = {{args.inopts}}
genome  = {{args.genome | quote}}
poscol  = 'genomic_pos_%s' % genome if genome!= 'hg38' else 'genomic_pos'

# get the genes
genes = genenorm(
	infile, 
	notfound = {{args.notfound | quote}},
	frm      = {{args.frm | quote}},
	to       = "symbol, genomic_pos_{{args.genome}}",
	genome   = genome,
	cachedir = {{args.cachedir | quote}},
	inopts   = inopts
)

reader = TsvReader(infile, **inopts)
writer = TsvWriter(outfile, ftype = 'bed')
for gene, hit in genes.items():
	
	# TODO: log those genes not found.
	if poscol not in hit: continue
	
	pos      = json.loads(hit[poscol])
	if not pos: continue 
	# TODO: have to figure out this (when a gene has isoforms)
	if isinstance(pos, list): pos = pos[0]
	
	r        = TsvRecord()
	r.CHR    = 'chr' + str(pos['chr'])
	if pos['strand'] == 1:
		r.STRAND = '+'
		r.START  = pos['start'] - {{args.up}}
		r.END    = pos['start'] + {{args.down}}
		r.END    = {% if args.incbody %}max(r.END, pos['end']){% endif %}
	else:
		r.STRAND = '-'
		r.END    = pos['end'] + {{args.up}}
		r.START  = pos['end'] - {{args.down}}
		r.START  = {% if args.incbody %}min(r.START, pos['start']){% endif %}
	
	r.NAME   = hit['symbol']
	r.SCORE  = '0'
	writer.write(r)
