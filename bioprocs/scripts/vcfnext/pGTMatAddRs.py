from bioprocs.utils import FileConn

{% python from bioprocs.utils import alwaysList %}
infile   = {{i.infile | quote}}
outfile  = {{o.outfile | quote}}
dbsnp    = {{args.dbsnp | quote}}
notfound = {{args.notfound | quote}}
exist    = {{args.exist | quote}}
chrorder = {{args.chrorder | alwaysList | repr}}
chrorder = dict(zip(chrorder, list(range(len(chrorder)))))

def compare(chr1, pos1, chr2, pos2):
	"""
	return 0 if they equal, -1 if snp1 less, else 1
	"""
	chr1 = chrorder.get(chr1, chr1)
	chr2 = chrorder.get(chr2, chr2)
	pos1 = int(pos1)
	pos2 = int(pos2)
	return -1 if (chr1, pos1) < (chr2, pos2) else 1 if (chr1, pos1) > (chr2, pos2) else 0

with FileConn(infile) as fin, FileConn(dbsnp) as fdb, FileConn(outfile, 'w') as fout:
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
				(chr_in, pos_in, rs, ref, alt) = names
			if exist == 'keep' and rs.startswith('rs'):
				fout.write(rin)
				rin = None
				continue
			(chr_db, pos_db, rs_db) = rdb.split('\t')[:3]
			c = compare(chr_in, pos_in, chr_db, pos_db)
			if c < 0:
				inparts[0] = '_'.join((chr_in, pos_in, notfound, ref, alt))
				fout.write('\t'.join(inparts))
				rin = None
				continue
			if c > 0:
				rdb = None
				continue
			inparts[0] = '_'.join((chr_in, pos_in, rs_db, ref, alt))
			fout.write('\t'.join(inparts))
			rin = None
		except StopIteration:
			break
		except:
			from traceback import print_exc
			print_exc()
			rin = rdb = None
			continue
