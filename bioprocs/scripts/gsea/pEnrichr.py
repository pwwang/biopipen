from os import path
from enrichr import Enrichr
from pyppl import Box
from bioprocs import rimport
from bioprocs.utils import shell2 as shell
from bioprocs.utils.parallel import Parallel
from bioprocs.utils.gene import genenorm
from bioprocs.utils.tsvio2 import TsvReader

{% from pyppl.utils import alwaysList %}
infile   = {{i.infile | quote}}
prefix   = {{i.infile | fn2 | quote}}
outdir   = {{o.outdir | quote}}
inopts   = {{args.inopts | repr}}
genecol  = {{args.genecol | ?bool | :_ | :0 | repr}}
top      = {{args.top}}
dbs      = {{args.libs | alwaysList | repr}}
plot     = {{args.plot | repr}}
nthread  = {{args.nthread | repr}}
Rscript  = {{args.Rscript | repr}}
cutoff   = {{args.cutoff | repr}}
devpars  = {{args.devpars | repr}}
pathview = {{args.pathview | repr}}
include  = {{args.include | ?bool | :_ | :None }}

shell.load_config(Rscript = Rscript)

if isinstance(cutoff, dict):
	if cutoff['by'] == 'p':
		cutoff['by'] = 'Pval'
	if cutoff['by'] == 'q':
		cutoff['by'] = 'AdjPval'

reader = TsvReader(infile, **inopts)
genes  = list(set([r[genecol] for r in reader if not include or include(r)]))
with open(path.join(outdir, 'genes.txt'), 'w') as fg:
	fg.write('\n'.join(genes))

en = Enrichr(cutoff = cutoff, top = top, Rscript = Rscript)
en.addList(genes, description = path.basename(infile))

para = Parallel(nthread = nthread)
runPathview = lambda r, hsa: shell.fg.Rscript(r, hsa)
for db in dbs:
	outfile = path.join(outdir, prefix + '.' + db + '.txt')
	en.enrich(db)
	en.export(outfile, top = 100)
	if plot:
		plotfile = path.join(outdir, prefix + '.' + db + '.png')
		en.plot(plotfile, res = devpars.res, width = devpars.width, height = devpars.height)
	if pathview and 'KEGG' in db:
		pathviewRDir  = path.join(outdir, prefix + '.' + db + '.pathview')
		pathviewRfile = path.join(pathviewRDir, 'pathview.R')
		shell.mkdir(pathviewRDir)
		with open(pathviewRfile, 'w') as f:
			f.write("""
			{rimport}('__init__.r')
			library(pathview)
			args = commandArgs(trailingOnly = TRUE)
			setwd({pathviewRDir!r})
			inopts = {{args.inopts | R}}
			inopts$rnames = FALSE
			indata = read.table.inopts({infile!r}, inopts)
			genes  = as.vector(indata[, {genecol}, drop = TRUE])
			pvargs = {{args.pathview | R}}
			{% raw %}
			if (!is.null(pvargs$fccol)) {{
				fcdata = as.vector(indata[, pvargs$fccol, drop = TRUE])
				names(fcdata) = genes
				genes = fcdata
			}}
			{% endraw %}
			pathview(gene.data = genes, pathway.id = args[1], species = 'hsa', gene.idtype="SYMBOL")
			""".format(
				rimport = rimport, genecol = genecol + 1 if isinstance(genecol, int) else genecol,
				infile = infile, pathviewRDir = pathviewRDir)
			)
		para.run(runPathview, [(pathviewRfile, term.Term.split('_')[-1]) for term in en.results[:top]])

