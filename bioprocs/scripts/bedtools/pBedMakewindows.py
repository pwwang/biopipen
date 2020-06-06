from diot import Diot
from bioprocs.utils import shell2 as shell

infile   = {{i.infile | quote}}
outfile  = {{o.outfile | quote}}
intype   = {{args.intype | quote}}
params   = {{args.params | repr}}
bedtools = {{args.bedtools | quote}}

shell.load_config(bedtools = bedtools)

if intype == 'bed':
	params.b = infile
else:
	params.g = infile

# params._out   = outfile
# params._debug = True
shell.bedtools.makewindows(**params).r > outfile
