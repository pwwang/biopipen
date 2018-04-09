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
	_outopts.update(outopts)
	outopts  = _outopts
	
	reader   = TsvReader(infile, **inopts)	
	if not reader.meta: reader.autoMeta()
	genecol  = genecol or reader.meta.keys()[0]
	
	genes    = list(set([r[genecol] for r in reader]))
	reader.rewind()
	
	dbfile   = path.join(cachedir, 'geneinfo.db')
	cache    = Cache(dbfile, 'geneinfo', {
		'_id': 'int primary key',
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
	
	wherelike = lambda x: ['|{}'.format(x), '|{}|%'.format(x), '%|{}'.format(x), '%|{}|%'.format(x)]
	foundfunc = lambda v, data: re.search(r'\|%s($|\|)' % v, data)
	factoryin = lambda v: ''.join(['|' + x for x in v]) if v and isinstance(v, list) else v if v else ''
	wherekeys = {
		'alias'          : ('alias[~~]', wherelike),
		'ensembl_protein': ('alias[~]', wherelike),
		'pfam'           : ('pfam[~]', wherelike),
		'uniprot'        : ('uniprot[~]', wherelike)
	}
	foundkeys = {
		'alias'          : foundfunc,
		'ensembl_protein': foundfunc,
		'pfam'           : foundfunc,
		'uniprot'        : foundfunc
	}
	cachefactory = {
		'alias'           : [factoryin] * 2,
		'ensembl_protein' : [factoryin] * 2,
		'pfam'            : [factoryin] * 2,
		'uniprot'         : [factoryin] * 2,
		'genomic_pos'     : [lambda x: json.dumps(x) if isinstance(x, (list, dict)) else str(x)] * 2,
		'genomic_pos_hg19': [lambda x: json.dumps(x) if isinstance(x, (list, dict)) else str(x)] * 2,
	}
	
	# query from cache
	tocols   = alwaysList(to)
	frmcols  = alwaysList(frm)
	columns  = list(set(tocols + frmcols))
	allfound = []
	allrest  = None
	for frm in alwaysList(frm):
		datafound, retrest = cache.query(columns, {frm: genes, 'taxid': [TAXIDS[genome]]*len(genes)}, wherekeys, foundkeys)
		allfound += datafound
		for key, val in retrest.items():
			if key == 'taxid': continue
			if allrest is None:
				allrest = set(val)
			else:
				allrest &= set(val)

	# query from api
	allrest = list(allrest)
	mgret = MyGeneInfo().querymany(allrest, scopes = frmcols, fields = columns, species = SPECIES[genome])
	
	# get all result for each query
	genetmp = {}
	for gene in mgret:
		if not gene['query'] in genetmp:
			genetmp[gene['query']] = []
		genetmp[gene['query']].append(gene)

	genemap = {}
	for query, gene in genetmp.items():
		
		# re-score the items if query is entirely matched
		for g in gene:
			# not all result returned
			if not all([x in g for x in tocols]): continue
			
			if any([g[x] == query for x in tocols]):
				g['_score'] += 10000
			elif any([str(g[x]).upper() == query.upper() for x in tocols]):
				g['_score'] += 5000
			elif any([x in g and query.upper() in [str(u).upper() for u in list(g[x]) for x in tocols]]):
				g['_score'] += 1000
		gene = sorted(
			[g for g in gene if '_score' in g], 
			key = lambda x: x['_score'], 
			reverse = True
		)
		if not gene: continue
		gene = gene[0]
		genemap[query] = {}
		for x in columns:
			if not x in gene:
				genemap[query][x] = None
			elif x in cachefactory \
				and isinstance(cachefactory[x], (tuple, list)) \
				and len(cachefactory[x]) > 1:
				genemap[query][x] = cachefactory[x][1](gene[x])
			else:
				genemap[query][x] = gene[x]
	del genetmp	

	# cache genemap
	cachedata = {}
	querys    = genemap.keys()
	for query in querys:
		for k, v in genemap[query].items():
			if not k in cachedata:
				cachedata[k] = []
			cachedata[k].append(v)
			
	if cachedata:
		cache.save(cachedata, cachefactory)
		del cachedata
	
	# restore cached data
	for gene in genes:
		# gene not cached
		if gene in allrest: continue
		# check if gene is the query
		datarow = None
		for row in allfound:
			for frmcol in frmcols:
				if frmcol in foundkeys and foundkeys[frmcol](gene, row[frmcol]):
					datarow = row
					break
				elif frmcol not in foundkeys and gene == row[frmcol]:
					datarow = row
					break
			if datarow: break
		if not datarow: continue
		genemap[gene] = datarow
	
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
			query = r[genecol]
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