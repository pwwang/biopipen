from pyppl import Box
from bioprocs.utils import logger
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter

{% python from os import path %}
{% python from pyppl.utils import alwaysList %}
{% assign cnames = dict(list = list, readlines = readlines,alwaysList = alwaysList, func = lambda x: x) %}
{% assign cntype = lambda x, path = path: 'none' if x is None else 'list' if isinstance(x, (list, tuple)) else 'readlines' if path.isfile(x) else 'func' if x.startswith('lambda') else 'alwaysList' %}
infile  = {{i.infile | quote}}
hfile   = {{i.hfile | quote}}
outfile = {{o.outfile | quote}}
inopts  = {{args.inopts | repr}}
cnames  = {{args.cnames | cnames.get(cntype(args.cnames), repr)}}

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
	reader.cnames = cnames(reader.cnames)

writer = TsvWriter(outfile, delimit = inopts.get('delimit', "\t"))
writer.cnames = reader.cnames
writer.writeHead()
writer.file.write(reader.file.read())

del reader
del writer
