import gzip

infile   = {{i.infile | quote}}
outfile  = {{o.outfile | quote}}
dbsnp    = {{args.dbsnp | quote}}
notfound = {{args.notfound | quote}}
chrsort  = {{args.chrsort | repr}}
if isinstance(chrsort, list):
	chrsort = dict(zip(chrsort, list(range(len(chrsort)))))

def mapChrForVersion(c):
	"""
	Only for human genome. If this is applied to other species, futher
	modification is needed.
	"""
	if c.startswith('chrM'):
		return 998
	elif c == 'chrX':
		return 999
	elif c == 'chrY':
		return 1000
	else:
		return int(c[3:])

def compare(chr1, pos1, chr2, pos2):
	"""
	return 0 if they equal, -1 if snp1 less, else 1
	"""
	pos1 = int(pos1)
	pos2 = int(pos2)
	if chrsort == 'version':
		chr1 = mapChrForVersion(chr1)
		chr2 = mapChrForVersion(chr2)
	elif chrsort == 'natural':
		pass # use original chr1, chr2
	else:
		chr1 = chrsort.get(chr1, chr1)
		chr2 = chrsort.get(chr2, chr2)
	return -1 if (chr1, pos1) < (chr2, pos2) else 1 if (chr1, pos1) > (chr2, pos2) else 0

with open(infile) as fin, \
	(gzip.open(dbsnp) if dbsnp.endswith('.gz') else open(dbsnp)) as fdb, \
	open(outfile, 'w') as fout:
	fout.write(fin.next())
	rin = rdb = None
	while True:
		try:
			rin = rin or fin.next()
			rdb = rdb or fdb.next()
			if rdb.startswith('#'):
				rdb = None
				continue
			inparts = rin.split('\t')
			names   = inparts[0].split('_')
			if len(names) == 4:
				(chr_in, pos_in, ref, alt) = names
			else:
				(chr_in, pos_in, _, ref, alt) = names
			(chr_db, pos_db, rs_db) = rdb.split('\t')[:3]
			cmp = compare(chr_in, pos_in, chr_db, pos_db)
			if cmp < 0:
				inparts[0] = '_'.join((chr_in, pos_in, notfound, ref, alt))
				fout.write('\t'.join(inparts))
				rin = None
				continue
			if cmp > 0:
				rdb = None
				continue
			inparts[0] = '_'.join((chr_in, pos_in, rs_db, ref, alt))
			fout.write('\t'.join(inparts))
			rin = rdb = None
		except StopIteration:
			break
		except:
			from traceback import print_exc
			print_exc()
			rin = rdb = None
			continue
