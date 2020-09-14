from os import path
from diot import Diot
from pyppl.utils import always_list
from bioprocs.utils import logger
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter

infile  = {{i.infile | quote}}
hfile   = {{i.hfile | quote}}
outfile = {{o.outfile | quote}}
inopts  = {{args.inopts | repr}}
cnames  = {{args.cnames | ?:isinstance(_, str) and _.startswith('lambda') | !repr}}

if isinstance(cnames, str):
	if path.isfile(cnames):
		with open(cnames) as f:
			cnames = [line.strip() for line in f.readlines() if line.strip()]
	else:
		cnames = always_list(cnames)
elif cnames and not callable(cnames):
	cnames = list(cnames)

if hfile and cnames and not callable(cnames):
	logger.warning('Both header file (hfile) and column names (args.cnames) specified, latter one discarded.')

header = None
if hfile:
	with open(hfile) as f:
		header = [h.strip() for h in f if h.strip()]

reader = TsvReader(infile, **inopts)
if header:
	reader.cnames = header
elif cnames and not callable(cnames):
	reader.cnames = cnames
if callable(cnames):
	reader.cnames = cnames(reader.cnames, path.basename(infile))

writer = TsvWriter(outfile, delimit = inopts.get('delimit', "\t"))
writer.cnames = reader.cnames
writer.writeHead()
writer.file.write(reader.file.read())

del reader
del writer
