from bioprocs.utils import runcmd, cmdargs

infile   = {{i.infile | repr}}
outfile  = {{o.outfile | repr}}
umfile   = {{o.umfile | repr}}
params   = {{args.params | repr}}
chain    = {{args.lochain | repr}}
liftover = {{args.liftover | repr}}

cmd = "{liftover} {infile} {chain} {outfile} {params}"
runcmd(cmd.format(
	liftover = liftover,
	infile   = str(repr(infile)),
	chain    = str(repr(chain)),
	outfile  = str(repr(outfile)),
	params   = cmdargs(params, dash = '-', equal = '=')
))
