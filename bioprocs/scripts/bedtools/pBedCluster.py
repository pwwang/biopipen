from diot import Diot
from bioprocs.utils import shell2 as shell

infile   = {{i.infile | quote}}
outfile  = {{o.outfile | quote}}
bedtools = {{args.bedtools | quote}}
params   = {{args.params | repr}}

shell.load_config(bedtools=bedtools)

params.i = infile
# params._stdout = outfile
shell.bedtools.cluster(**params).r > outfile
