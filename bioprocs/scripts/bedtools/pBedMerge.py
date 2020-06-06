from diot import Diot
from bioprocs.utils import shell2 as shell

params = {{args.params | repr}}
params.i = {{i.infile | quote}}

shell.load_config(bedtools = {{args.bedtools | quote}})

shell.bedtools.merge(**params).r > {{out.outfile | squote}}
