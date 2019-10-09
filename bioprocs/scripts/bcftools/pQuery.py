from pyppl import Box
from bioprocs.utils import shell2 as shell

infile   = {{i.infile | quote}}
outfile  = {{o.outfile | quote}}
bcftools = {{args.bcftools | quote}}
params   = {{args.params | repr}}

shell.load_config(bcftools = bcftools)

params._ = infile
params._out = outfile
shell.bcftools.query(**params)
