
infile  = {{i.infile | repr}}
outfile = {{o.outfile | repr}}

chrom = {{args.chr | quote}}

with open(infile) as f, open(outfile, 'w') as fout:
	for line in f:
		if line.startswith('#'):
			fout.write(line)
		elif line.startswith('chr'):
			fout.write(line)
		else:
			fout.write(chrom + line)
