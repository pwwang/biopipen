from collections import defaultdict
from pyppl import Box
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter

infile  = {{ i.infile | quote }}
outfile = {{ o.outfile | quote }}
inopts  = {{args.inopts | repr}}
infmt   = {{args.infmt | quote}}
cutoff  = {{args.cutoff | repr}}

degrees = defaultdict(lambda: 0)
if infmt.startswith('pair'):
	reader = TsvReader(infile, **inopts)
	for r in reader:
		if cutoff:
			try:
				score = float(r[2])
			except TypeError:
				raise TypeError('The 3rd column should be a score for apply the cutoff.')
			if score < cutoff:
				continue
		degrees[r[0]] += 1
		degrees[r[1]] += 1
	writer = TsvWriter(outfile)
	for node in sorted(degrees.keys(), key = lambda x: degrees[x], reverse = True):
		if infmt.endswith('complete'):
			writer.write([node, int(int(degrees[node])/2)])
		else:
			writer.write([node, degrees[node]])
	writer.close()
else:
	raise ValueError('Input format other than "pair" not supported yet.')