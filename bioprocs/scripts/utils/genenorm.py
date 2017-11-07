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
	def genenorm(infile, outfile = None, col = 0, notfound = 'ignore', frm = 'symbol, alias', to = 'symbol', header = False, genome = 'human', skip = 0, delimit = '\t', tmpdir = gettempdir(), comment = '#'):
		species = {
			'hg19': 'human',
			'hg38': 'human',
			'mm9' : 'mouse',
			'mm10': 'mouse'
		}
		species  = species[genome] if genome in species else genome
		cachedir = path.join(tmpdir, 'mygene-cache')
		if not path.isdir(cachedir): makedirs(cachedir)

		skip += int(header)
		with open(infile) as f:
			genes   = sorted(set([line.split(delimit)[col] for i, line in enumerate(f.read().splitlines()) if i >= skip and line.strip() and not line.startswith(comment)]))
		
		uid   = md5(''.join(genes) + str(frm) + str(to) + str(species)).hexdigest()
		cache = path.join(cachedir, 'mygeneinfo.%s' % uid)
		if path.isfile(cache):
			with open(cache) as f: mgret = json.load(f)
		else:
			mgret = MyGeneInfo().getgenes(genes, scopes = frm, fields = to, species = species, size=1)
			with open(cache, 'w') as f:	json.dump(mgret, f, indent = 4)
		
		genemap = {}
		for gene in mgret:
			if to not in gene: continue
			genemap[gene['query']] = gene[to]

		if not outfile: return genemap

		with open(infile) as f, open(outfile, 'w') as fout:
			i = 0
			for line in f:
				if i < skip:
					fout.write(line)
					i += 1
					continue

				line2 = line.strip()
				if not line2 or line2.startswith(comment):
					fout.write(line)
					continue

				parts = line2.split(delimit)
				gene  = parts[col]
				if not gene in genemap:
					if notfound == 'error':
						raise ValueError('Cannot convert %s to %s.' % (gene, to))
					if notfound == 'skip':
						continue
					if notfound == 'ignore':
						genemap[gene] = gene
				parts[col] = genemap[gene]
				fout.write(delimit.join(parts) + '\n')
		return genemap, cache