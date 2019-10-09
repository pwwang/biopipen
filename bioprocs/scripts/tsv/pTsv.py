from pyppl import Box
from bioprocs.utils.tsvio2 import TsvWriter, TsvReader, TsvRecord

infile  = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
inopts  = {{args.inopts | repr}}
outopts = {{args.outopts | repr}}
helper  = {{args.helper | repr}}
if not isinstance(helper, list):
	helper = [helper]

helper = [line for line in helper if line]
exec('\n'.join(helper), globals())

row = {{args.row}}
inopts['row'] = {{args.inopts.get('row', None)}}

reader       = TsvReader(infile, **inopts)
writer       = TsvWriter(outfile, delimit = outopts.get('delimit', "\t"))
outcnames    = outopts.get('cnames', True)
headCallback = {{args.outopts.get('headCallback', True)}}
if outcnames is True:
	writer.cnames = reader.cnames
	writer.writeHead(headCallback)
elif isinstance(outcnames, (list, tuple)):
	writer.cnames = list(outcnames)
	writer.writeHead(headCallback)

for r in reader:
	rec = row(r)
	if rec is False:
		continue
	if rec is None or rec is True:
		rec = r
	writer.write(rec)
writer.close()

