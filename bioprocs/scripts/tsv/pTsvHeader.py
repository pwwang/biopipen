from pyppl import Box
from bioprocs.utils.tsvio2 import TsvReader

infile  = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
inopts  = {{args.inopts | repr}}
cfilter = {{args.filter}}

reader = TsvReader(infile, **inopts)
cnames = reader.cnames
if callable(cfilter):
	cnames = cfilter(cnames)
with open(outfile, 'w') as f:
	f.write("".join([cname + "\n" for cname in cnames if cname]))
