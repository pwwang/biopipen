from os import path
from pyppl import Box
from bioprocs.utils import shell

infile   = {{i.infile | quote}}
outfile  = {{o.outfile | quote}}
params   = {{args.params | repr}}
bedtools = {{args.bedtools | quote}}
gsize    = {{args.gsize | quote}}

if not gsize:
	raise ValueError('Genome size file is required (args.gsize).')
if not path.isfile(gsize):
	raise ValueError('Genome size file does not exist (args.gsize).')

bedtools = shell.Shell({'bedtools': bedtools}, subcmd = True, equal = ' ', dash = '-').bedtools

params.i = infile
params.g = gsize
params._stdout = outfile
bedtools.shift(**params).run()