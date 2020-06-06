from diot import Diot
from bioprocs.utils import shell2 as shell

infile   = {{i.infile | quote}}
outfile  = {{o.outfile | quote}}
params   = {{args.params | repr}}
bedtools = {{args.bedtools | quote}}

shell.load_config(bedtools= bedtools)

params.ibam    = infile
# params._stdout = outfile
shell.bedtools.genomecov(**params).r > outfile
