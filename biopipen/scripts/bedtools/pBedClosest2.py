from diot import Diot
from bioprocs.utils import shell2 as shell

afile    = {{i.afile | quote}}
bfiles   = {{i.bfiles | repr}}
outfile  = {{o.outfile | quote}}
params   = {{args.params | repr}}
bedtools = {{args.bedtools | quote}}

shell.load_config(bedtools= bedtools)

params.a = afile
params.b = bfiles
# params._stdout = outfile
shell.bedtools.closest(**params).r > outfile
