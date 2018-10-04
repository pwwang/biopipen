from pyppl import Box
from bioprocs.utils.gene import genenorm
from bioprocs.utils.tsvio import TsvReader, TsvWriter

infile   = {{i.infile | quote}}
outfile  = {{o.outfile | quote}}
notfound = {{args.notfound | quote}}
inopts   = {{args.inopts}}
outopts  = {{args.outopts}}
genecol  = {{args.genecol | repr}}
fromipi  = {{args.fromipi | bool}}
ipidb    = {{args.ipidb | quote}}

if genecol is None:
	genecol = 0

# get dicts 
ipireader = TsvReader(ipidb, ftype = 'nometa')
dicts = {}
for r in ipireader:
	gene = r[6]
	ipis = r[9].split(';')
	if not gene: continue
	if fromipi:
		dicts.update({ipi:gene for ipi in ipis if ipi and gene})
	else:
		dicts[gene] = '|'.join(ipi for ipi in ipis if ipi)

reader = TsvReader(infile, **inopts)
writer = TsvWriter(outfile, **outopts)
writer.meta = reader.meta
if outopts.query:
	writer.meta.add('_QUERY')

if outopts.head:
	writer.writeHead(
		delimit   = outopts.headDelimit,
		prefix    = outopts.headPrefix,
		transform = outopts.headTransform)

for r in reader:
	query = r[genecol]
	if not query in dicts:
		if notfound == 'ignore':
			ret = query
		elif notfound == 'skip':
			continue
		elif notfound == 'error':
			raise ValueError('{} not found!'.format(query))
	else:
		ret = dicts[query]
	r[genecol] = ret
	if outopts.query:
		if isinstance(r, list):
			r.append(query)
		else:
			r._QUERY = query
	writer.write(r)
	
