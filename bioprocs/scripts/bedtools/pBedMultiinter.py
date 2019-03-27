from pyppl import Box
from bioprocs.utils import shell

infiles  = {{i.infiles | repr}}
outfile  = {{o.outfile | quote}}
params   = {{args.params | repr}}
bedtools = {{args.bedtools | quote}}

bedtools = shell.Shell({'bedtools': bedtools}, subcmd = True, equal = ' ', dash = '-').bedtools

params.i       = infiles
params._stdout = outfile
bedtools.multiinter(**params).run()
