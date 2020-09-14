from diot import Diot
from biopipen.utils import shell2 as shell

afile    = {{i.afile | quote}}
bfile    = {{i.bfile | quote}}
outfile  = {{o.outfile | quote}}
params   = {{args.params | repr}}
bedtools = {{args.bedtools | quote}}

shell.load_config(bedtools= bedtools)


params.a = afile
params.b = bfile
# params._stdout = outfile
shell.bedtools.closest(**params).r > outfile
