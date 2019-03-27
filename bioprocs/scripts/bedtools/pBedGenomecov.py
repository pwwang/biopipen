from pyppl import Box
from bioprocs.utils import shell

infile   = {{i.infile | quote}}
outfile  = {{o.outfile | quote}}
params   = {{args.params | repr}}
bedtools = {{args.bedtools | quote}}

bedtools = shell.Shell({'bedtools': bedtools}, subcmd = True, equal = ' ', dash = '-').bedtools

params.ibam    = infile
params._stdout = outfile
bedtools.genomecov(**params).run()