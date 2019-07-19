from pyppl import Box
from bioprocs.utils import shell2 as shell

infile  = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
params  = {{args.params | repr}}
inopts  = {{args.inopts | repr}}
case    = {{args.case | repr}}
mem     = {{args.mem | repr}}
tmpdir  = {{args.tmpdir | quote}}
unique  = {{args.unique | repr}}

shell.load_config(
	sort = dict(_dupkey = True, _env = {'LANG': 'C'} if case else {'LANG': 'en_US.UTF-8'}))

params.t = params.get('t', inopts.get('delimit', "\t"))
params.T = tmpdir
params.u = unique
params.S = mem

if 'cnames' in inopts:
	inopts.skip += 1

if inopts.skip:
	shell.head(_ = infile, n = inopts.skip, _out = outfile)
	#params.__stdout = outfile
	shell.pipe.tail(_ = infile, n = '+' + str(inopts.skip + 1)) | \
		shell.sort(_out_ = outfile, **params)
else:
	params._    = infile
	params._out = outfile
	shell.sort(**params)
