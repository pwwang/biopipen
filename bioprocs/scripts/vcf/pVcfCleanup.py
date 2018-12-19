
from os import path
from bioprocs.utils.tsvio import TsvReader

infile  = {{i.infile | repr}}
outfile = {{o.outfile | repr}}
ref     = {{args.ref | repr}}
refdict = {{args.ref | prefix | quote}} + '.dict'
reffai  = {{args.ref | quote}} + '.fai'

if not path.isfile(refdict) and not path.isfile(reffai):
	raise OSError('A dict file or fai file not exists for the reference file.')

contigs = []
if path.isfile(refdict):
	reader = TsvReader(refdict, skip = 1)
	reader.meta.add('FLAG', 'CONFIG')
	for r in reader:
		contigs.append(r.CONFIG[3:])

elif path.isfile(reffai):
	reader = TsvReader(reffai)
	reader.meta.add('CONFIG')
	for r in reader:
		contigs.append(r.CONFIG)

with open(infile) as f, open(outfile, 'w') as fout:
	for line in f:
		if line.startswith('#'):
			fout.write(line)
		else:
			parts = line.split('\t')
			if parts[0] in contigs:
				fout.write(line)
