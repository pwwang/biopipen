from os import path
from diot import Diot
from bioprocs.utils import shell2 as shell

infile   = {{i.infile | quote}}
outfile  = {{o.outfile | quote}}
params   = {{args.params | repr}}
bedtools = {{args.bedtools | quote}}
gsize    = {{args.gsize | quote}}

if not gsize:
	raise ValueError('Genome size file is required (args.gsize).')
if not path.isfile(gsize):
	raise ValueError('Genome size file does not exist (args.gsize).')

shell.load_config(bedtools=bedtools)

params.i = infile
params.g = gsize
# params._stdout = outfile
shell.bedtools.shift(**params).r > outfile
