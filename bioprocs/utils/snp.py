from os import path
from sys import stderr
from cruzdb import Genome
from tempfile import gettempdir
from pyppl import Box
from bioprocs.utils import alwaysList
from bioprocs.utils.cache import Cache
from bioprocs.utils.tsvio import TsvReader, TsvWriter, TsvRecord

class RecordNotFound(Exception):
	pass

def snpinfo(infile, outfile = None, notfound = 'ignore', genome = 'hg19', dbsnpver = "150", inopts = None, outopts = None, snpcol = None, cachedir = gettempdir()):
	_inopts = Box(skip = 0, comment = '#', delimit = '\t')
	_inopts.update(inopts or {})
	inopts  = _inopts

	_outopts = Box(delimit = '\t', headDelimit = '\t', headPrefix = '', headTransform = None, head = True, ftype = 'bed', cnames = 'refUCSC, alleles, alleleFreqs, alleleFreqCount')
	_outopts.update(outopts)
	outopts  = _outopts
	cnames   = alwaysList(outopts['cnames'])

	reader = TsvReader(infile, **inopts)
	if not reader.meta: reader.autoMeta()
	snpcol = snpcol or reader.meta.keys()[0]

	snps = list(set([r[snpcol] for r in reader]))
	reader.rewind()

	dbfile = path.join(cachedir, 'snpinfo_%s_%s.db' % (genome, dbsnpver))
	schema = {
		'chrom'          : 'text',              # chr8
		'chromStart'     : 'int',               # 128700232L
		'chromEnd'       : 'end',               # 128700233L
		'name'           : 'text primary key',  # rs7005394
		'score'          : 'real',              # 0
		'strand'         : 'text',              # +
		'refNCBI'        : 'text',              # T
		'refUCSC'        : 'text',              # T
		'observed'       : 'text',              # C/T
		'class'          : 'single',            # single
		'avHet'          : 'real',              # 0.49
		'avHetSE'        : 'real',              # 0.02
		'func'           : 'text',              # set(['ncRNA'])
		'submitterCount' : 'int',               # 20
		'submitters'     : 'text',              # 1000GENOMES,ABI,BCM-HGSC-SUB,BCM_SSAHASNP,BGI,BL ...
		'alleleFreqCount': 'int',               # 2
		'alleles'        : 'text',              # C,T,
		'alleleNs'       : 'text',              # 2634.000000,2374.000000,
		'alleleFreqs'    : 'text',              # 0.525958,0.474042,
	}
	cache   = Cache(dbfile, 'snpinfo', schema, 'name')
	dummies = {
		'func': dict(
			query  = Cache.DUMMY['array']['query'],
			find   = Cache.DUMMY['array']['find'],
			insert = lambda col, data: (col, ' // ' + ' // '.join(Cache._uniqueData(list(data), True))),
			update = lambda col, data: (col, Function.concat(Field(col), value = ' // ' + ' // '.join(Cache._uniqueData(list(d), True)))),
			result = lambda data: data if isinstance(data, list) else list(filter(None, data.split(' // ')))
		)
	}
	columns = [
		'chrom', 'name', 'chromStart', 'chromEnd', 'score', 'strand'
	] + cnames
	ret, allrest = cache.query(columns, {'name': snps}, dummies)

	dbsnp = Genome(db = genome)
	dbsnp = getattr(dbsnp, "snp%s" % dbsnpver)

	writer = None
	if outfile:
		head          = outopts['head']
		headPrefix    = outopts['headPrefix']
		headDelimit   = outopts['headDelimit']
		headTransform = outopts['headTransform']
		del outopts['head']
		del outopts['headPrefix']
		del outopts['headDelimit']
		del outopts['headTransform']
		writer = TsvWriter(outfile, **outopts)
		if head:
			writer.writeHead(prefix = headPrefix, delimit = headDelimit, transform = headTransform)

	if writer:
		for r in ret.values():
			r.CHR    = r.chrom
			r.START  = r.chromStart
			r.END    = r.chromEnd
			r.NAME   = r.name
			r.SCORE  = r.score
			r.STRAND = r.strand
			writer.write(r)

	cached = []
	if allrest:
		for snp in allrest['name']:
			s = dbsnp.filter_by(name=snp).first()
			if not s:
				if notfound == 'error':
					raise RecordNotFound('Record not found: %s' % snp)
				elif notfound == 'skip':
					continue
				else:
					stderr.write('Record not found: %s \n' % snp)
					continue
			cached.append(s)
			if writer:
				r = TsvRecord()
				r.CHR    = s.chrom
				r.START  = s.chromStart
				r.END    = s.chromEnd
				r.NAME   = s.name
				r.SCORE  = s.score
				r.STRAND = s.strand
				for cname in cnames:
					setattr(r, cname, getattr(s, cname))
				writer.write(r)

	# save cached data
	cachedata = {}
	for c in cached:
		for k in schema.keys():
			if not k in cachedata:
				cachedata[k] = []
			cachedata[k].append(getattr(c, k))
	if cachedata:
		cache.save(cachedata, dummies)
	return {r.name: r for r in ret.values() + cached}
