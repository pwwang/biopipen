from pyppl import Box
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter
infile  = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
ncol    = {{args.ncol}}
nameopt = {{args.name | quote}}
inopts  = {{args.inopts | repr}}

if ncol not in [3,6,8,33,66,88]:
	raise ValueError('Unknown ncol, expect one of [3,6,8,66,88]')

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
		return [chrom, int(pos) - 1, pos]
	if ncol == 6:
		return [chrom, int(pos) - 1, pos, name, 0, '+']
	if ncol == 8:
		return [chrom, int(pos) - 1, pos, name, 0, '+', ref, alt]
	if ncol == 33:
		return [chrom, int(pos) - 1, pos] + row[1:]
	if ncol == 66:
		return [chrom, int(pos) - 1, pos, name, 0, '+'] + row[1:]
	if ncol == 88:
		return [chrom, int(pos) - 1, pos, name, 0, '+', ref, alt] + row[1:]

dft_inopts = Box(cnames = True, attach = False, row = rowFactory)
dft_inopts.update(inopts)

reader = TsvReader(infile, **dft_inopts)
writer = TsvWriter(outfile)
for r in reader:
	writer.write(r)
writer.close()


