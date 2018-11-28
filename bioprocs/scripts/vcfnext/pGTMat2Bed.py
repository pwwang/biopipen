from bioprocs.utils.tsvio2 import TsvReader, TsvWriter
infile  = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
ncol    = {{args.ncol}}
nameopt = {{args.name | quote}}

def rowFactory(row):
	rparts = row[0].split('_')
	if len(rparts) == 4:
		(chrom, pos, ref, alt) = rparts
		name = chrom + '_' + pos if nameopt == 'neat' else row[0]
	elif len(rparts) == 5:
		(chrom, pos, name, ref, alt) = rparts
		if name == 'NOVEL':
			name = chrom + '_' + pos
		if nameopt == 'full':
			name = row[0]
	else:
		raise ValueError('Malformat genotype matrix, expect 4 or 5 items in row names.')
	if ncol == 3:
		return [chrom, pos, int(pos) + 1]
	if ncol == 6:
		return [chrom, pos, int(pos) + 1, name, 0, '+']
	else:
		return [chrom, pos, int(pos) + 1, name, 0, '+', ref + ',' + alt, ','.join(rows[1:])]

reader = TsvReader(infile, cnames = True, attach = False, row = rowFactory)
writer = TsvWriter(outfile)
for r in reader:
	writer.write(r)
writer.close()


