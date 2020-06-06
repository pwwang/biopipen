from diot import Diot
from bioprocs.utils import shell2 as shell

infile   = {{i.infile | repr}}
outfile  = {{o.outfile | repr}}
umfile   = {{o.umfile | repr}}
params   = {{args.params | str}}
chain    = {{args.lochain | repr}}
liftover = {{args.liftover | repr}}

shell.load_config(liftover = liftover)

shell.liftover(infile, chain, outfile, umfile, **params).fg
