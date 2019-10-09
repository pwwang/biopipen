from pyppl import Box
from bioprocs.utils import funcargs
from bioprocs.utils.tsvio2 import TsvWriter, TsvRecord
from gff import Gff

infile    = {{i.infile | quote}}
outfile   = {{o.outfile | quote}}
bedcols   = {{args.bedcols | repr}}
keepattrs = {{args.keepattrs | repr}}
outhead   = {{args.outhead | repr}}

bedcols.NAME = bedcols.get('NAME', 'lambda attrs: \
										attrs["id"] if "id" in attrs else \
										attrs["name"] if "name" in attrs else \
										attrs["CHR"] + ":" + attrs["START"] + "-" + attrs["END"]')
# convert strings to functions
for key in bedcols:
	bedcols[key] = eval(bedcols[key])


writer = TsvWriter(outfile)
writer.cnames = ['CHR', 'START', 'END', 'NAME', 'SCORE', 'STRAND']
writer.cnames.extend([field for field in bedcols if field != 'NAME'])

if keepattrs:
	writer.cnames.append('ATTRIBUTES')
	bedcols.ATTRIBUTES = lambda attrs: ';'.join(
		['{0}={1}'.format(key, val) for key, val in attrs.items()
		 if key not in writer.cnames])

if outhead:
	writer.writeHead(lambda cns: ('' if outhead is True else outhead) + '\t'.join(cns))

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
	for key in bedcols:
		r[key] = bedcols[key](attrs)

	writer.write(r)
writer.close()
