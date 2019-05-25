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

writer = TsvWriter(outfile)
if inopts.cnames:
	writer.cnames = reader.cnames + ['Pathway_' + lib for lib in libs]
	writer.writeHead()

reader.rewind()
for r in reader:
	gene = r[genecol]
	for lib in libs:
		if not lib in pathways:
			pathways[lib] = {}
		if gene not in pathways[lib]:
			en = Enrichr(lib)
			gmap = en.genemap(gene)
			pathways[lib][gene] = '|'.join(gmap.terms if gmap else [])
		r['Pathway_' + lib] = pathways[lib][gene]
	writer.write(r.values())

reader.close()
writer.close()
