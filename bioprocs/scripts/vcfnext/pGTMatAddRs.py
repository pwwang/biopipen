import gzip

infile   = {{i.infile | quote}}
outfile  = {{o.outfile | quote}}
dbsnp    = {{args.dbsnp | quote}}
notfound = {{args.notfound | quote}}

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
			inparts                     = rin.split('\t')
			(chr_in, pos_in, ref, alt)  = inparts[0].split('_')
			(chr_db, pos_db, rs_db)     = rdb.split('\t')[:3]
			if (chr_in, int(pos_in)) < (chr_db, int(pos_db)):
				inparts[0] = '_'.join((chr_in, pos_in, notfound, ref, alt))
				fout.write('\t'.join(inparts))
				rin = None
				continue
			if (chr_in, int(pos_in)) > (chr_db, int(pos_db)):
				rdb = None
				continue
			inparts[0] = '_'.join((chr_in, pos_in, rs_db, ref, alt))
			fout.write('\t'.join(inparts))
			rin = rdb = None
		except StopIteration:
			break
		except:
			continue
