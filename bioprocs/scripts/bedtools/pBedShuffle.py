from pyppl import Box
from bioprocs.utils import shell

infile   = {{ i.infile | quote}}
outfile  = {{ o.outfile | quote}}
n        = {{ args.n | repr}}
gsize    = {{ args.gsize | quote}}
seed     = {{ args.seed | repr}}
params   = {{ args.params | repr}}
bedtools = {{ args.bedtools | quote}}
bedtools = shell.Shell({'bedtools': bedtools}, subcmd = True, dash = '-', equal = ' ').bedtools

if n and n<1:
	total = shell.wcl(infile)
	n = int(total * n)

params.i = infile
params.g = gsize
params.seed = seed

if n:
	shuffle = bedtools.shuffle(**params)
	shuffle.pipe().head(n = n, _stdout = outfile)
else:
	params._stdout = outfile
	shuffle = bedtools.shuffle(**params)
	shuffle.run()
