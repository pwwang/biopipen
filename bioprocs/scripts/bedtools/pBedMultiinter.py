from diot import Diot
from bioprocs.utils import shell2 as shell

infiles  = {{i.infiles | repr}}
outfile  = {{o.outfile | quote}}
params   = {{args.params | repr}}
bedtools = {{args.bedtools | quote}}
shell.load_config(bedtools= bedtools)

params.i       = infiles
# params._stdout = outfile
shell.bedtools.multiinter(**params).r > outfile
