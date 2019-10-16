from os import path
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter

infile  = {{i.infile | quote}}
samfile = {{i.samfile | quote}}
outfile = {{o.outfile | quote}}
samples = {{args.samples | repr}}

if samfile and path.isfile(samfile):
	with open(samfile) as f:
		samples = f.read().strip().splitlines()
elif samfile:
	samples = samfile.split(',')
elif not isinstance(samples, list):
	samples = samples.split(',')

reader = TsvReader(infile)
writer = TsvWriter(outfile)
writer.cnames = reader.cnames
writer.writeHead()

for r in reader:
	if r.Tumor_Sample_Barcode in samples:
		writer.write(r)

reader.close()
writer.close()
