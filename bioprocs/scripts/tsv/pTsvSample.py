from pyppl import Box
from bioprocs.utils import shell

infile    = {{i.infile | quote}}
outfile   = {{o.outfile | quote}}
inopts    = {{args.inopts | repr}}
n         = {{args.n | repr}}
arsample  = {{args.arsample | quote}}
replace   = {{args.replace | repr}}
keeporder = {{args.keeporder |repr}}
seed      = {{args.seed |repr}}
params    = {{args.params | repr}}
arsample  = shell.Shell(dict(sample = arsample)).sample

if inopts.get('skip', 0):
	shell.head(n = inopts.skip, _ = infile, _stdout = outfile)
	infile_skipped = outfile + '.skipped'
	shell.tail(n = '+' + str(inopts.skip + 1), _ = infile, _stdout = infile_skipped)
	infile = infile_skipped

params._ = infile
params._stdout_ = outfile
params.k = n
params.d = seed
if keeporder:
	params.p = True
if replace:
	params.r = True
else:
	params.o = True

arsample(**params).run()
