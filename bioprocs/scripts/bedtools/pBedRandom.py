from pyppl import Box
from bioprocs.utils import shell

l        = {{i.l | repr}}
n        = {{i.n | repr}}
seed     = {{args.seed | repr}}
bedtools = {{args.bedtools | quote}}
gsize    = {{args.gsize | quote}}
outfile  = {{o.outfile | quote}}

shell.TOOLS.bedtools = bedtools
shell.Shell(subcmd = True, dash = '-', equal = '=').bedtools.random(l = l, n = n, seed = seed, g = gsize, _stdout = outfile)
