from diot import Diot
from bioprocs.utils import shell2 as shell

l        = {{i.l | repr}}
n        = {{i.n | repr}}
seed     = {{args.seed | repr}}
bedtools = {{args.bedtools | quote}}
gsize    = {{args.gsize | quote}}
outfile  = {{o.outfile | quote}}

shell.load_config(bedtools = bedtools)
shell.bedtools.random(l = l, n = n, seed = seed, g = gsize).r > outfile
