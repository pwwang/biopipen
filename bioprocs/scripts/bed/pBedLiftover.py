from bioprocs.utils import runcmd, cmdargs

infile   = {{in.infile | repr}}
outfile  = {{out.outfile | repr}}
umfile   = {{out.umfile | repr}}
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
