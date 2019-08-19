import re, json, requests
from os import path
from mygene import MyGeneInfo
from medoo import Raw, Field
from pyppl import Box
from pyppl.utils import alwaysList
from bioprocs.utils.cache import Cache
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter, TsvRecord
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
# local to remote
FIELD_L2M = {
	'ensembl_gene'      : 'ensembl.gene',
	'ensembl_protein'   : 'ensembl.protein',
	'ensembl_transcript': 'ensembl.transcript',
	'refseq_genomic'    : 'refseq.genomic',
	'refseq_rna'        : 'refseq.rna',
	'refseq_protein'    : 'refseq.protein',
	'uniprot_Swiss_Prot': 'uniprot.Swiss-Prot',
}
# remote to local
FIELD_M2L = {
	'ensembl.gene'      : 'ensembl_gene',
	'ensembl.protein'   : 'ensembl_protein',
	'ensembl.transcript': 'ensembl_transcript',
	'refseq.genomic'    : 'refseq_genomic',
	'refseq.rna'        : 'refseq_rna',
	'refseq.protein'    : 'refseq_protein',
	'uniprot.Swiss-Prot': 'uniprot_Swiss_Prot'
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
	rets  = []
	try:
		mgret = MyGeneInfo().querymany(*args, **kwargs)
	except requests.exceptions.ConnectionError:
		return rets
	for ret in mgret:
		out = {}
		rets.append(out)
		for key, val in ret.items():
			if 'ensembl' == key:
				ensembl = val[0] if isinstance(val, list) else (val or {})
				out['ensembl_gene']       = ensembl.get('gene', '')
				out['ensembl_protein']    = ensembl.get('protein', [])
				out['ensembl_transcript'] = ensembl.get('transcript', [])
			elif 'refseq' == key:
				refseq = val[0] if isinstance(val, list) else (val or {})
				out['refseq_genomic'] = refseq.get('genomic', [])
				out['refseq_rna']     = refseq.get('rna', [])
				out['refseq_protein'] = refseq.get('protein', [])
			elif 'uniprot' == key:
				uniprot = val[0] if isinstance(val, list) else (val or {})
				out['uniprot_Swiss_Prot'] = uniprot.get('Swiss-Prot', [])
			else:
				out[key] = val
	return rets

fields2local  = lambda keys: [FIELD_M2L.get(key, key) for key in keys]
fields2remote = lambda keys: [FIELD_L2M.get(key, key) for key in keys]

def genenorm(infile, outfile = None, notfound = 'ignore', frm = 'symbol, alias', to = 'symbol', genome = 'hg19', inopts = None, outopts = None, genecol = None, cachedir = gettempdir()):

	_inopts = Box(skip = 0, comment = '#', delimit = '\t')
	_inopts.update(inopts or {})
	inopts  = _inopts

	_outopts = Box(delimit = '\t', append = False, query = False, head = True)
	_outopts.update(outopts or {})
	outopts  = _outopts
	outquery = outopts.get('query', False)
	outhead  = outopts.get('head', outopts.get('cnames', True))
	if 'query' in outopts:
		outquery = outopts['query']
		del outopts['query']
	if 'head' in outopts:
		outhead = outopts['head']
		del outopts['head']
	if 'cnames' in outopts:
		outhead = outopts['cnames']
		del outopts['cnames']

	reader  = TsvReader(infile, **inopts)
	#if not reader.meta: reader.autoMeta()
	genecol = genecol or 0
	genes   = set()
	ncol    = 0
	for r in reader:
		ncol = ncol or len(r)
		genes.add(r[genecol].strip())
	reader.rewind()
	if not reader.meta:
		reader.meta.extend(['COL' + str(i + 1) for i in range(ncol)])
	genes = list(genes)

	dbfile   = path.join(cachedir, 'geneinfo.db')
	cache    = Cache(dbfile, 'geneinfo', {
		'_id'               : 'text primary key',
		'symbol'            : 'text',
		'HGNC'              : 'int',
		'alias'             : "text default ''",
		'ensembl_gene'      : 'text',
		'ensembl_protein'   : 'text',
		'ensembl_transcript': 'text',
		'refseq_genomic'    : 'text',
		'refseq_rna'        : 'text',
		'refseq_protein'    : 'text',
		'entrezgene'        : 'int',
		'genomic_pos'       : 'text',
		'genomic_pos_hg19'  : 'text',
		'genomic_pos_mm9'   : 'text',
		'ipi'               : 'text',
		'pfam'              : "text default ''",
		'pdb'               : 'text',
		'type_of_gene'      : 'text',
		'taxid'             : 'int',
		'uniprot_Swiss_Prot': "text default ''",
	}, '_id')

	dummies = {
		'symbol'            : 'iplain',
		'alias'             : 'iarray',
		'pfam'              : 'iarray',
		'uniprot'           : 'iarray',
		'genomic_pos'       : 'json',
		'genomic_pos_hg19'  : 'json',
		'genomic_pos_mm9'   : 'json',
		'ipi'               : 'iarray',
		'pdb'               : 'iarray',
		'refseq_genomic'    : 'iarray',
		'refseq_protein'    : 'iarray',
		'refseq_rna'        : 'iarray',
		'ensembl_protein'   : 'iarray',
		'ensembl_transcript': 'iarray',
		'uniprot_Swiss_Prot': 'iarray',
	}

	# query from cache
	tocols = alwaysList(to)
	# alias
	tocols = replaceList(tocols, ['ensg', 'ensemblgene', 'ensembl'], 'ensembl.gene')
	tocols = replaceList(tocols, ['uniprot'], 'uniprot.Swiss-Prot')
	tocols = replaceList(tocols, ['refseq'], 'refseq.rna')
	tocols = fields2local(tocols)

	frmcols = alwaysList(frm)
	frmcols = replaceList(frmcols, ['ensg', 'ensemblgene', 'ensembl'], 'ensembl.gene')
	frmcols = replaceList(frmcols, ['uniprot'], 'uniprot.Swiss-Prot')
	frmcols = replaceList(frmcols, ['refseq'], 'refseq.rna')
	frmcols = fields2local(frmcols)

	columns  = list(set(tocols + frmcols + ['taxid']))
	frmkeys  = ','.join(frmcols)
	allfound, allrest = cache.query(columns, {frmkeys: genes, 'taxid': TAXIDS[genome]}, dummies)

	# query from api
	mgret = querygene(allrest[frmkeys], scopes = fields2remote(frmcols), fields = fields2remote(columns), species = SPECIES[genome])
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
		writer = TsvWriter(outfile, **outopts)
		writer.meta.extend(reader.meta)
		if outquery:
			write.meta.append('_QUERY')
			
		if len(tocols) > 1:
			gcolidx = genecol if isinstance(genecol, int) else writer.meta.index(genecol)
			writer.meta[(gcolidx+1):(gcolidx+1)] = [(tocol, None) for tocol in tocols[1:]]
		
		if outhead:
			writer.writeHead()
		#print writer.meta

		#i = 0
		for row in reader:
			r = TsvRecord(row.values(), reader.meta)

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
				#if (i <= 10): print genecol
				r[genecol] = genemap[query][tocols[0]]
				#if (i <= 10): print genemap[query][tocols[0]], r
				if len(tocols) > 1:
					for tocol in tocols[1:]:
						r[tocol] = genemap[query][tocol]

			if outquery:
				r._QUERY = query

			#if (i <= 10): print r
			i += 1
			writer.write(r)

	return genemap
