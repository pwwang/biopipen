from diot import Diot
from bioprocs.utils import shell2 as shell

infile   = {{ i.infile | quote}}
outfile  = {{ o.outfile | quote}}
n        = {{ args.n | repr}}
gsize    = {{ args.gsize | quote}}
seed     = {{ args.seed | repr}}
params   = {{ args.params | repr}}
bedtools = {{ args.bedtools | quote}}

shell.load_config(bedtools=bedtools)

if n and n<1:
	total = shell.wcl(infile)
	n = int(total * n)

params.i = infile
params.g = gsize
params.seed = seed

if n:
	shell.bedtools.shuffle(**params).p | shell.head(n=n) ^ outfile
else:
	shell.bedtools.shuffle(**params).r > outfile
