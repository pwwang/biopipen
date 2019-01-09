from bioprocs.utils.tsvio2 import TsvWriter, TsvRecord
from gff import Gff

infile    = {{i.infile | quote}}
outfile   = {{o.outfile | quote}}
attr2name = {{args.attr2name}}
keepinfo  = {{args.keepinfo | repr}}

writer = TsvWriter(outfile)
writer.cnames = ['CHR', 'START', 'END', 'NAME', 'SCORE', 'STRAND']
if keepinfo:
	writer.cnames.append('ORIGINAL')

def getNameFromAttrs(attrs):
	if attr2name:
		return attr2name(**attrs)
	for key in sorted(attrs.keys()):
		if key in writer.cnames:
			continue
		if 'id' in key.lower():
			return attrs[key]
		if 'name' in key.lower():
			return attrs[key]
		return attrs[key]

gff = Gff(infile)
for record in gff:
	r        = TsvRecord()
	r.CHR    = record['seqid']
	r.START  = record['start']
	r.END    = record['end']
	r.SCORE  = record['score']
	r.STRAND = record['strand']
	attrs    = record['attributes']
	attrs.update(dict(
		CHR    = r.CHR,
		START  = r.START,
		END    = r.END,
		SCORE  = r.SCORE,
		STRAND = r.STRAND
	))
	r.NAME   = getNameFromAttrs(attrs)
	if keepinfo:
		r.ORIGINAL = '; '.join('{}={}'.format(k,v) for k, v in attrs.items() if k not in writer.cnames)
	writer.write(r)
	
