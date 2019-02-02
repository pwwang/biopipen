from bioprocs.utils import shell

infile   = {{i.infile | repr}}
outfile  = {{o.outfile | repr}}
umfile   = {{o.umfile | repr}}
params   = {{args.params | repr}}
chain    = {{args.lochain | repr}}
liftover = {{args.liftover | repr}}

shell.TOOLS.liftover = liftover
shell.Shell(dash = '-', equal = '=').liftover(infile, chain, outfile, **params).run()

