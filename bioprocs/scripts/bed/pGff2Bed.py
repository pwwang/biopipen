from bioprocs.utils.tsvio import TsvWriter, TsvRecord
from gff import Gff

infile    = {{in.infile | quote}}
outfile   = {{out.outfile | quote}}
attr2name = {{args.attr2name | lambda a: 'None' if a is None else a}}
keepinfo  = {{args.keepinfo | repr}}

writer = TsvWriter(outfile, ftype = 'bed')
if keepinfo:
	writer.meta.add('SOURCE', 'TYPE', 'PHASE', 'ATTRS')

def getNameFromAttrs(attrs):
	if attr2name:
		return attr2name(**attrs)
	else:
		for key, val in attrs.items():
			if 'id' in key.lower():
				return val
			elif 'name' in key.lower():
				return val
			else:
				return val

gff = Gff(infile)
for record in gff:
	r        = TsvRecord()
	r.CHR    = record['seqid']
	r.START  = record['start']
	r.END    = record['end']
	r.SCORE  = record['score']
	r.STRAND = record['strand']
	attrs    = record['attributes']
	r.NAME   = getNameFromAttrs(attrs)
	if keepinfo:
		r.SOURCE = record['source']
		r.TYPE   = record['type']
		r.PHASE  = record['phase']
		r.ATTRS  = ';'.join(['{k}={v}'.format(k=k, v=v) for k, v in attrs.items()])
	writer.write(r)
	
