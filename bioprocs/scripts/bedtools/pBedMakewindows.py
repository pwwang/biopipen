from pyppl import Box
from bioprocs.utils import shell

infile    = {{i.infile | quote}}
outfile  = {{o.outfile | quote}}
intype    = {{args.intype | quote}}
params   = {{args.params | repr}}
bedtools = {{args.bedtools | quote}}

bedtools = shell.Shell({'bedtools': bedtools}, subcmd = True, equal = ' ', dash = '-').bedtools

if intype == 'bed':
	params.b = infile
else:
	params.g = infile

params._       = intfile
params._stdout = outfile
bedtools.makewindows(**params).run()