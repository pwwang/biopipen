from os import path
from enrichr import Enrichr
from pyppl import Box
from bioprocs.utils.gene import genenorm
from bioprocs.utils.tsvio2 import TsvReader

{% python from bioprocs.utils import alwaysList %}
infile   = {{i.infile | quote}}
prefix   = {{i.infile | fn2 | quote}}
outdir   = {{o.outdir | quote}}
inopts   = {{args.inopts}}
genecol  = {{args.genecol or 0 | repr}}
top      = {{args.top}}
dbs      = {{args.libs | alwaysList | repr}}
plot     = {{args.plot | repr}}
devpars  = {{args.devpars | repr}}
cutoff   = {{args.cutoff | repr}}
if isinstance(cutoff, dict):
	if cutoff['by'] == 'p':
		cutoff['by'] = 'Pval'
	if cutoff['by'] == 'q':
		cutoff['by'] = 'AdjPval'

reader = TsvReader(infile, **inopts)
genes  = [r[genecol] for r in reader]

en = Enrichr(cutoff = cutoff, top = top)
en.addList(genes, description = path.basename(infile))

for db in dbs:
	outfile = path.join(outdir, prefix + '.' + db + '.txt')
	en.enrich(db)
	en.export(outfile, top = 100)
	if plot:
		plotfile = path.join(outdir, prefix + '.' + db + '.png')
		en.plot(plotfile, res = devpars.res, width = devpars.width, height = devpars.height)
