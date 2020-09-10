from diot import Diot
from bioprocs.utils import shell2 as shell
from bioprocs.utils.reference import vcfIndex

infiles  = {{i.infiles | repr}}
outfile  = {{o.outfile | quote}}
nthread  = {{args.nthread | repr}}
bcftools = {{args.bcftools | quote}}
params   = {{args.params | repr}}
gz       = {{args.gz | repr}}
tabix    = {{args.tabix | repr}}

shell.load_config(bcftools = bcftools)

infiles = [vcfIndex(infile, tabix) for infile in infiles]

params._       = infiles
params.o       = outfile
params.threads = nthread
params.O       = 'z' if gz else 'v'
shell.bcftools.concat(**params).fg()
