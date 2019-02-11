from pyppl import Box
from bioprocs.utils import shell
from bioprocs.utils.reference import vcfIndex

vcffile = {{ i.vcffile | quote}}
regfile = {{ i.regfile | quote}}
outfile = {{ o.outfile | quote}}
tabix   = {{ args.tabix | quote}}
params  = {{ args.params | repr}}

shell.TOOLS.tabix = tabix
vcffile = vcfIndex(vcffile, tabix)

params._ = [vcffile, regfile]
params._stdout = outfile
shell.Shell().tabix(**params).run()

