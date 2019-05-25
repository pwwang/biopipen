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
enrich   = Enrichr()
pathways = {}

writer = TsvWriter(outfile)
if inopts.cnames:
	writer.cnames = reader.cnames + ['Pathway_' + lib for lib in libs]
	writer.writeHead()

reader.rewind()
for r in reader:
	gene = r[genecol]
	if gene not in pathways:
		pathways[gene] = {lib:enlib.terms
			for lib, enlib in enrich.genemap(gene, libs).items()}

	for lib in libs:
		r['Pathway_' + lib] = ' | '.join(pathways[gene][lib]) \
			if lib in pathways[gene] else ''
	writer.write(r.values())

reader.close()
writer.close()
