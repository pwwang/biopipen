import re, json
from os import path
from mygene import MyGeneInfo
from medoo import Raw, Field
from pyppl import Box
from bioprocs.utils import alwaysList
from bioprocs.utils.cache import Cache
from bioprocs.utils.tsvio import TsvReader, TsvWriter
from tempfile import gettempdir

"""
`notfound`: What if a symbol is not found. Default: skip
	- skip  : skip the record(don't write it to output file)
	- ignore: use the original name;
	- error : report erro
"""
SPECIES = {
	'hg19': 'human',
	'hg38': 'human',
	'mm9' : 'mouse',
	'mm10': 'mouse'
}
TAXIDS  = {
	'hg19': 9606,
	'hg38': 9606,
	'mm9' : 10090,
	'mm10': 10090
}
class RecordNotFound(Exception):
	pass

def genenorm(infile, outfile = None, notfound = 'ignore', frm = 'symbol, alias', to = 'symbol', genome = 'hg19', inopts = None, outopts = None, genecol = None, cachedir = gettempdir()):

	_inopts = Box(skip = 0, comment = '#', delimit = '\t')
	_inopts.update(inopts or {})
	inopts  = _inopts

	_outopts = Box(delimit = '\t', headDelimit = '\t', headPrefix = '', headTransform = None, head = True, query = False)
	_outopts.update(outopts or {})
	outopts  = _outopts

	reader   = TsvReader(infile, **inopts)
	if not reader.meta: reader.autoMeta()
	genecol  = genecol or reader.meta.keys()[0]
	genes    = list(set([r[genecol].strip() for r in reader]))
	reader.rewind()

	dbfile   = path.join(cachedir, 'geneinfo.db')
	cache    = Cache(dbfile, 'geneinfo', {
		'_id': 'text primary key',
		'symbol': 'text',
		'HGNC': 'int',
		'alias': 'text',
		'ensembl': 'text',
		'ensembl_protein': 'text',
		'entrezgene': 'int',
		'genomic_pos': 'text',
		'genomic_pos_hg19': 'text',
		'pfam': 'text',
		'taxid': 'int',
		'uniprot': 'text'
	}, '_id')

	dummies = {
		'symbol' : 'iplain',
		'alias'  : 'iarray',
		'pfam'   : 'iarray',
		'uniprot': 'iarray',
		'genomic_pos': 'json',
		'genomic_pos_hg19': 'json',
	}

	# query from cache
	tocols   = alwaysList(to)
	frmcols  = alwaysList(frm)
	frmkeys  = ','.join(frmcols)
	columns  = list(set(tocols + frmcols + ['taxid']))
	allfound, allrest = cache.query(columns, {frmkeys: genes, 'taxid': TAXIDS[genome]}, dummies)

	# query from api
	mgret = MyGeneInfo().querymany(allrest[frmkeys], scopes = frmcols, fields = columns, species = SPECIES[genome])
	# get all result for each query
	genetmp = {}
	for gret in mgret:
		if not gret['query'] in genetmp:
			genetmp[gret['query']] = []
		genetmp[gret['query']].append(gret)

	genemap   = {}
	data2save = {}
	for query, gret in genetmp.items():
		# re-score the items if query is entirely matched
		score = 0
		gr    = None
		for g in gret:
			# not all result returned
			if not all([x in g for x in tocols]): continue

			if any([g[x] == query for x in tocols]):
				thescore = g['_score'] + 10000
			elif any([str(g[x]).upper() == query.upper() for x in tocols]):
				thescore = g['_score'] + 5000
			elif any([x in g and query.upper() in [str(u).upper() for u in list(g[x]) for x in tocols]]):
				thescore = g['_score'] + 1000
			else:
				thescore = g['_score']
			if thescore > score:
				score = thescore
				gr    = g

		if not gr: continue
		del gr['_score']
		del gr['query']

		gr = Cache._result({x:(gr[x] if x in gr else '') for x in set(columns + gr.keys())}, dummies)
		for x, val in gr.items():
			if not x in data2save:
				data2save[x] = []
			data2save[x].append(val)
		genemap[query] = gr

	# add cached data
	for i, ret in allfound.items():
		query = genes[i]
		genemap[query] = ret
	#del genetmp
	#print genemap

	# cache genemap
	#cachedata = {}
	#querys    = genemap.keys()
	#for query in querys:
	#	for k, v in genemap[query].items():
	#		if not k in cachedata:
	#			cachedata[k] = []
	#		cachedata[k].append(v)

	#if cachedata:
	#	cache.save(cachedata, cachefactory)
	#	del cachedata
	if data2save:
		cache.save(data2save, dummies)

	if outfile:
		writer   = TsvWriter(outfile, delimit = outopts['delimit'])
		writer.meta.update(reader.meta)
		if len(tocols) > 1:
			items   = writer.meta.items()
			gcolidx = writer.meta.keys().index(genecol)
			items[(gcolidx+1):(gcolidx+1)] = [(tocol, None) for tocol in tocols[1:]]
			writer.meta.clear()
			writer.meta.add(*items)
		if outopts['query']:
			writer.meta.add('_QUERY')

		if outopts['head']:
			writer.writeHead(outopts['headPrefix'], outopts['headDelimit'], outopts['headTransform'])
		for r in reader:
			query = r[genecol].strip()
			if query not in genemap:
				if notfound == 'error':
					raise RecordNotFound('Record not found: %s' % query)
				elif notfound == 'skip':
					continue
				if len(tocols) > 1:
					for tocol in tocols[1:]:
						r[tocol] = ''
			else:
				r[genecol] = genemap[query][tocols[0]]
				if len(tocols) > 1:
					for tocol in tocols[1:]:
						r[tocol] = genemap[query][tocol]

			if outopts['query']:
				r._QUERY = query
			writer.write(r)

	return genemap
