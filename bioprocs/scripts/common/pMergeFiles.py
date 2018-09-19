
infiles = {{i.infiles | repr}}
outfile = {{o.outfile | quote}}
header  = {{args.header | repr}}

with open(infiles.pop(0)) as fin, open(outfile, 'w') as fout:
	fout.write(fin.read())

	for infile in infiles:
		with open(infile) as f:
			if header:
				f.readline()
			fout.write(f.read())
			