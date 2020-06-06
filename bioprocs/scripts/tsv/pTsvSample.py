from diot import Diot
from bioprocs.utils import shell2 as shell

infile    = {{i.infile | quote}}
outfile   = {{o.outfile | quote}}
inopts    = {{args.inopts | repr}}
n         = {{args.n | repr}}
arsample  = {{args.arsample | quote}}
replace   = {{args.replace | repr}}
keeporder = {{args.keeporder |repr}}
seed      = {{args.seed | repr}}
params    = {{args.params | repr}}

shell.load_config(arsample = arsample)

if inopts.get('skip', 0):
	shell.head(n = inopts.skip, _ = infile).r > outfile
	infile_skipped = outfile + '.skipped'
	shell.tail(n = '+' + str(inopts.skip + 1), _ = infile).r > infile_skipped
	infile = infile_skipped

params._    = infile
params._out = outfile
params.k    = n
if seed:
	params.d    = seed
if keeporder:
	params.p = True
if replace:
	params.r = True
else:
	params.o = True

shell.arsample(**params).fg
