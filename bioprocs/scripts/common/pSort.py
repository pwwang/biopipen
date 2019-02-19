from pyppl import Box
from bioprocs.utils import shell

infile  = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
params  = {{args.params | repr}}
inopts  = {{args.inopts | repr}}
case    = {{args.case | repr}}
mem     = {{args.mem | repr}}
tmpdir  = {{args.tmpdir | quote}}
unique  = {{args.unique | repr}}

shellargs = dict(dash = '-', equal = ' ', duplistkey = True)
shellargs['env'] = {'LANG': 'C'} if case else {'LANG': 'en_US.UTF-8'}

params.t = params.get('t', inopts.get('delimit', "\t"))
params.T = tmpdir
params.u = unique
params.S = mem

if inopts.skip:
	shell.head(_ = infile, n = inopts.skip, _stdout = outfile)
	#params.__stdout = outfile
	shell.Shell().tail(_ = infile, n = '+' + str(inopts.skip + 1)) \
		.pipe(**shellargs).sort(__stdout = outfile, **params).run()
else:
	params._ = infile
	params._stdout = outfile
	shell.Shell(**shellargs).sort(**params).run()
