from diot import Diot
from bioprocs.utils import shell

infile   = {{i.infile | quote}}
outfile  = {{o.outfile | quote}}
bedtools = {{args.bedtools | quote}}
params   = {{args.params | repr}}

bedtools = shell.Shell(dict(bedtools = bedtools), subcmd = True).bedtools

params.i = infile
params._stdout = outfile
bedtools(**params).run()
