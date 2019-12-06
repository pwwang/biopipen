
from bioprocs.utils import gztype, shell2 as shell
from bioprocs.utils.reference import vcfIndex

infile  = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
outidx  = {{o.outidx | quote}}
tabix   = {{args.tabix | quote}}

ofile = vcfIndex(infile, tabix)
shell.ln_s(ofile, outfile)
shell.ln_s(ofile + '.tbi', outidx)
