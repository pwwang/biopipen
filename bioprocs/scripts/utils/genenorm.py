# PYPPL REPEAT START: genenorm
if 'genenorm' not in vars() or not callable(genenorm):
	from os import path, makedirs
	from tempfile import gettempdir
	from hashlib import md5
	from mygene import MyGeneInfo
	import json
	"""
	`notfound`: What if a symbol is not found. Default: skip
		- skip  : skip the record(don't write it to output file)
		- ignore: use the original name;
		- error : report erro
	"""
	def genenorm(infile, outfile = None, notfound = 'ignore', frm = 'symbol, alias', to = 'symbol', genome = 'human', inopts = None, outopts = None, inmeta = ['GENE'], genecol = None, tmpdir = gettempdir()):
		species = {
			'hg19': 'human',
			'hg38': 'human',
			'mm9' : 'mouse',
			'mm10': 'mouse'
		}
		inopts1  = inopts or {}
		inopts   = {'skip': 0, 'comment': '#', 'delimit': '\t'}
		inopts.update(inopts1)
		outopts1 = outopts or {}
		outopts  = {'headprefix': '#', 'delimit': '\t', 'meta': False, 'head': True, 'metaprefix': '##META/'}
		outopts.update(outopts1)

		species  = species[genome] if genome in species else genome
		cachedir = path.join(tmpdir, 'mygene-cache')
		if not path.isdir(cachedir): makedirs(cachedir)

		genes   = []
		if isinstance(inmeta, list):
			greader = readBase(infile, **inopts)
			greader.meta.add(*inmeta)
		else:
			greader = globals()['read' + inmeta[0].upper() + inmeta[1:]](infile, **inopts)

		genecol = genecol or 'GENE'

		for r in greader:
			g = r[genecol]
			if not g in genes:
				genes.append(g)

		genes      = sorted(genes)
		alwaysList = lambda x: x if isinstance(x, list) else [y.strip() for y in x.split(',')]
		frm        = alwaysList(frm)
		to         = alwaysList(to)
		to2        = sorted(list(set(frm) | set(to)))
		uid        = md5(''.join(genes) + str(frm) + str(to) + str(species)).hexdigest()
		cache      = path.join(cachedir, 'mygeneinfo.%s' % uid)

		if path.isfile(cache):
			with open(cache) as f: mgret = json.load(f)
		else:
			mgret = MyGeneInfo().querymany(genes, scopes = frm, fields = to2, species = species)
			with open(cache, 'w') as f:	json.dump(mgret, f, indent = 4)
		
		genetmp = {}
		for gene in mgret:
			if not gene['query'] in genetmp:
				genetmp[gene['query']] = []
			genetmp[gene['query']].append(gene)
		
		genemap = {}
		for query, gene in genetmp.items():
			# re-score the items if query is entirely matched
			for g in gene:
				if not all([x in g for x in to]): continue
				
				if any([g[x] == query for x in to]):
					g['_score'] += 10000
				elif any([str(g[x]).upper() == query.upper() for x in to]):
					g['_score'] += 5000
				elif any([x in g and query.upper() in [str(u).upper() for u in list(g[x]) for x in to]]):
					g['_score'] += 1000
			gene = sorted([g for g in gene if any([x in g for x in to])], key = lambda x: x['_score'], reverse = True)
			if not gene: continue
			if len(to) == 1:
				genemap[query] = gene[0][to[0]]
			else:
				genemap[query] = {x: (gene[0][x] if x in gene[0] else None) for x in to}

		del genetmp

		if not outfile: return genemap, None
		greader.rewind()
		gwriter = writeBase(outfile)
		gwriter.meta.borrow(greader.meta)
		if len(to) > 0:
			gwriter.meta.add(*to)

		if outopts['meta']:
			gwriter.writeMeta(prefix = outopts['metaprefix'])
		if outopts['head']:
			gwriter.writeHead(prefix = outopts['headprefix'], delimit = outopts['delimit'])

		for r in greader:
			g = r[genecol]
			for t in to: 
				if t == genecol: continue
				r[t] = ''

			if g not in genemap:
				if notfound == 'error':
					raise ValueError('Cannot convert %s to %s.' % (g, to))
				elif notfound == 'skip':
					continue
				#else: # ignore, don't change r.GENE
				#	pass
			else:
				newg = genemap[g]
				if isinstance(newg, dict):
					for k,v in newg.items():
						r[k] = v
					r[genecol] = newg[to[0]]
				else:
					r[to[0]] = newg
					r[genecol] = newg
			gwriter.write(r, delimit = outopts['delimit'])

		return genemap, cache
# PYPPL REPEAT END: genenorm
