from os import path
from pyppl import Box
from bioprocs.utils import shell2 as shell

infile  = {{i.infile | quote}}
outdir  = {{o.outdir | quote}}
jvarkit = {{args.jvarkit | quote}}
params  = {{args.params | repr}}
Rscript = {{args.Rscript | quote}}

shell.load_config(jvarkit = jvarkit, Rscript = Rscript)

plotR = path.join(outdir, 'plot.R')
with open(plotR, 'w') as f:
	f.write('setwd(%r)\n' % outdir)

params._ = infile
params._raise = True
params._out = plotR
shell.jvarkit.vcfstatsjfx(**params)


shell.fg.Rscript(plotR)
