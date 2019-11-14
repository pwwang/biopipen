
from os import path
from bioprocs.utils import FileConn
from bioprocs.utils.tsvio2 import TsvReader

infile  = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
ref     = {{args.ref | quote}}
reffai  = ref + '.fai'
refdict = ref[:-3] + '.dict'

if not path.isfile(refdict) and not path.isfile(reffai):
	raise OSError('A dict file or fai file not exists for the reference file.')

contigs = []
if path.isfile(refdict):
	reader = TsvReader(refdict, skip = 1, cnames = False)
	for r in reader:
		contigs.append(r[1][3:])

elif path.isfile(reffai):
	reader = TsvReader(reffai, cnames = False)
	for r in reader:
		contigs.append(r[0])

with FileConn(infile, 'rt') as f, open(outfile, 'w') as fout:
	for line in f:
		# also remove contigs from header
		if line.startswith('##contig=<ID='):
			contig = line[13:].split(',')[0].split('>')[0]
			if contig in contigs:
				fout.write(line)
		elif line.startswith('#'):
			fout.write(line)
		else:
			config, _ = line.split('\t', 1)
			if config in contigs:
				fout.write(line)
