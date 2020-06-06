from diot import Diot
from bioprocs.utils import shell2 as shell
from bioprocs.utils.reference import vcfIndex

vcffile = {{ i.vcffile | quote}}
regfile = {{ i.regfile | quote}}
outfile = {{ o.outfile | quote}}
tabix   = {{ args.tabix | quote}}
params  = {{ args.params | repr}}

shell.load_config(tabix=tabix)
vcffile = vcfIndex(vcffile, tabix)

params._ = [vcffile, regfile]
# params._stdout = outfile
shell.tabix(**params).r > outfile
