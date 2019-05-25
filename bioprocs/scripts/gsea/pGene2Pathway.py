from pyppl import Box
from enrichr import Enrichr
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter

{%python from bioprocs.utils import alwaysList%}
infile  = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
inopts  = {{args.inopts | repr}}
genecol = {{args.genecol | repr}}
libs    = {{args.libs | alwaysList}}

reader   = TsvReader(infile, **inopts)
genes    = reader.dump(genecol)
genes    = list(set(genes))
pathways = {}
for lib in libs:
	en = Enrichr(lib)
	pathways[lib] = {}
	for gene in genes:
		gmap = en.genemap(gene)
		pathways[lib][gene] = '|'.join(gmap.terms if gmap else [])

writer = TsvWriter(outfile)
if inopts.cnames:
	writer.cnames = reader.cnames + ['Pathway_' + lib for lib in libs]
	writer.writeHead()

reader.rewind()
for r in reader:
	for lib in libs:
		r['Pathway_' + lib] = pathways[lib][r[genecol]]
	writer.write(r.values())

reader.close()
writer.close()
