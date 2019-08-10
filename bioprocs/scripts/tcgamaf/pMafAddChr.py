
from pyppl import Box
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter

infile  = {{i.infile | quote}}
outfile = {{o.outfile | quote}}

# Keep the comments
with open(infile) as fin, open(outfile, 'w') as fout:
	for line in fin:
		if line.startswith('#'):
			fout.write(line)
		else:
			break

reader = TsvReader(infile)
writer = TsvWriter(outfile, append = True)
writer.cnames = reader.cnames
writer.writeHead()
for r in reader:
	if not r.Chromosome.startswith('chr'):
		r.Chromosome = 'chr' + r.Chromosome
	writer.write(r)
reader.close()
writer.close()
