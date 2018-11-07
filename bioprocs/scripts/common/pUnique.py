from pyppl import Box
from bioprocs.utils.tsvio import TsvReader, TsvWriter

infile  = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
inopts  = {{args.inopts}}
col     = {{args.col | repr}}
outopts  = {
	'head': False, 
	'headPrefix': '', 
	'headDelimit': '\t', 
	'headTransform': None, 
	'delimit': '\t',
	'ftype': '',
	'cnames': ''
}
outopts.update({{args.outopts}})

head = outopts['head']
headPrefix = outopts['headPrefix']
headDelimit = outopts['headDelimit']
headTransform = outopts['headTransform']
del outopts['head']
del outopts['headPrefix']
del outopts['headDelimit']
del outopts['headTransform']

reader  = TsvReader(infile, **inopts)
if not reader.meta:
	reader.autoMeta()
writer  = TsvWriter(outfile, **outopts)
if not writer.meta:
	writer.meta.update(reader.meta)
if head:
	writer.writeHead(delimit = headDelimit, prefix = headPrefix, transform = headTransform)

colkeys = []
if col != '*':
	if isinstance(col, int):
		colkeys = [reader.meta.keys()[col]]
	elif isinstance(col, (tuple, list)):
		colkeys = [reader.meta.keys()[c] for c in col]
		
# try to keep the original order
{% if args.sorted %}
last = None
for r in reader:
	identity = str(r) if not colkeys else str([r[c] for c in colkeys])
	if last is None or last != identity:
		last = identity
		writer.write(r)
{% else %}
pool = {}
for r in reader:
	identity = str(r) if not colkeys else str([r[c] for c in colkeys])
	if not identity in pool:
		pool[identity] = 1
		writer.write(r)
{% endif %}