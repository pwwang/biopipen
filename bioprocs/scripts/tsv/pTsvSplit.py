from diot import Diot
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter

infile  = {{i.infile | quote}}
outdir  = {{o.outdir | quote}}
prefix  = {{i.infile | stem2 | quote}}
ext     = {{i.infile | ext | quote}}
inopts  = {{args.inopts | repr}}
outopts = {{args.outopts | repr}}
by      = {{args.by | ? :isinstance(_, str) and _.startswith('col:')
					| = [4:]
     				| $repr }}

reader    = TsvReader(infile, **inopts)
outheader = outopts.pop('header', True)

by_type = None
if isinstance(by, int):
	by_type = 'size'
elif isinstance(by, str):
	by_type = 'column'
	if by.isdigit():
		by = int(by)
elif callable(by):
	by_type = 'callable'
else:
	raise ValueError('Unsupported `args.by = %r`' % by)

def get_writer(tag):
	outfile = outdir + "/" + prefix + "_" + str(tag) + ext
	writer = TsvWriter(outfile, **outopts)
	if outheader is True:
		writer.cnames = reader.cnames
		if writer.cnames:
			writer.writeHead()
	elif outheader:
		writer.cnames = list(outheader)
		writer.writeHead()
	return writer

outfiles = {}
for i, r in enumerate(reader):
	if by_type == 'size':
		tag = int(i / by)
	elif by_type == 'column':
		tag = r[by]
	else:
		tag = by(r)
	if tag not in outfiles:
		outfiles[tag] = get_writer(tag)
	outfiles[tag].write(r)

reader.close()
for writer in outfiles.values():
	writer.close()
