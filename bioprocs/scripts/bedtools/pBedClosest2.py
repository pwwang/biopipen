from pyppl import Box
from bioprocs.utils import shell

afile    = {{i.afile | quote}}
bfiles   = {{i.bfiles | repr}}
outfile  = {{o.outfile | quote}}
params   = {{args.params | repr}}
bedtools = {{args.bedtools | quote}}

bedtools = shell.Shell({'bedtools': bedtools}, subcmd = True, equal = ' ', dash = '-').bedtools

params.a = afile
params.b = bfiles
params._stdout = outfile
bedtools.closest(**params).run()