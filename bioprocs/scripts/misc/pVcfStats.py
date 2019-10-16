from pyppl import Box
from bioprocs.utils import shell2 as shell

infile  = {{i.infile | quote}}
outdir  = {{o.outdir | quote}}
jvarkit = {{args.jvarkit | quote}}
params  = {{args.params | repr}}

shell.load_config(jvarkit = jvarkit)

params._ = infile
params.o = outdir
params._raise = True
shell.fg.jvarkit.vcfstats(**params)

shell.sed(i = 's/rotate by 90 right/"rotate by 90 right"/', _ = outdir + '/Makefile')
shell.make(_cwd = outdir, _raise = True)
