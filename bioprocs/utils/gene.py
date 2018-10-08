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

def replaceList(l, search, replace):
	if not isinstance(search, list):
		search = [search]
	ret = l[:]
	for i, e in enumerate(ret):
		if e in search:
			ret[i] = replace
	return ret

def querygene(*args, **kwargs):
	# change ensg to ensemblgene
	kwargs['scopes'] = replaceList(kwargs['scopes'], ['ensg', 'ensembl.gene', 'ensembl'], 'ensemblgene')
	kwargs['fields'] = replaceList(kwargs['fields'], ['ensg', 'ensemblgene', 'ensembl'], 'ensembl.gene')
	
	if 'ensemblgene' in kwargs['scopes'] and not 'ensembl.gene' in kwargs['fields']:
		kwargs['fields'].append('ensembl.gene')
	mgret = MyGeneInfo().querymany(*args, **kwargs)
	for ret in mgret:
		if 'ensembl' in ret:
			if 'gene' in ret['ensembl']:
				ret['ensembl'] = ret['ensembl']['gene']
			# more than one ensembl gene mapped, take the first one
			elif isinstance(ret['ensembl'], list) and len(ret['ensembl']) > 0 and 'gene' in ret['ensembl'][0]:
				ret['ensembl'] = ret['ensembl'][0]['gene']
	return mgret

def genenorm(infile, outfile = None, notfound = 'ignore', frm = 'symbol, alias', to = 'symbol', genome = 'hg19', inopts = None, outopts = None, genecol = None, cachedir = gettempdir()):

	_inopts = Box(skip = 0, comment = '#', delimit = '\t')
	_inopts.update(inopts or {})
	inopts  = _inopts

	_outopts = Box(delimit = '\t', headDelimit = '\t', headPrefix = '', headTransform = None, head = True, query = False)
	_outopts.update(outopts or {})
	outopts  = _outopts

	reader   = TsvReader(infile, **inopts)
	if not reader.meta: reader.autoMeta()
	genecol  = genecol or 0
	genes    = list(set([r[genecol].strip() for r in reader]))
	reader.rewind()

	dbfile   = path.join(cachedir, 'geneinfo.db')
	cache    = Cache(dbfile, 'geneinfo', {
		'_id': 'text primary key',
		'symbol': 'text',
		'HGNC': 'int',
		'alias': "text default ''",
		'ensembl': 'text',
		'ensembl_protein': 'text',
		'entrezgene': 'int',
		'genomic_pos': 'text',
		'genomic_pos_hg19': 'text',
		'pfam': "text default ''",
		'taxid': 'int',
		'uniprot': "text default ''",
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
	tocols   = replaceList(tocols, ['ensg', 'ensemblgene', 'ensembl.gene'], 'ensembl')
	frmcols  = alwaysList(frm)
	frmcols  = replaceList(frmcols, ['ensg', 'ensemblgene', 'ensembl.gene'], 'ensembl')
	frmkeys  = ','.join(frmcols)
	columns  = list(set(tocols + frmcols + ['taxid']))
	allfound, allrest = cache.query(columns, {frmkeys: genes, 'taxid': TAXIDS[genome]}, dummies)
	# query from api
	mgret = querygene(allrest[frmkeys], scopes = frmcols, fields = columns, species = SPECIES[genome])
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
		# make it unique
		ds_keys = list(data2save.keys())
		data2save_uniq = {k:[] for k in ds_keys}
		tmp_container  = []
		for i in range(len(data2save[ds_keys[0]])):
			tmp = {k:data2save[k][i] for k in ds_keys}
			if not tmp in tmp_container:
				tmp_container.append(tmp)
				for k in ds_keys:
					data2save_uniq[k].append(data2save[k][i])

		cache.save(data2save_uniq, dummies)

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

		#i = 0
		for r in reader:
			i += 1

			#if (i <= 10): print r
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
				if isinstance(r, list):
					r.append(query)
				else:
					r._QUERY = query

			#if (i <= 10): print r
			writer.write(r)

	return genemap
