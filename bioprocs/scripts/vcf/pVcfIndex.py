from os import path
from bioprocs.utils import gztype, shell2 as shell
from bioprocs.utils.reference import vcfIndex

infile  = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
outidx  = {{o.outidx | quote}}
tabix   = {{args.tabix | quote}}

def to_outdir(filepath, out):
	if path.islink(filepath):
		shell.ln_s(filepath, out)
	else:
		# if it is not a link, it will be cleaned, as it is in <indir>
		# we move it to <outdir> to keep it available for cache check
		shell.mv(filepath, out)

ofile = vcfIndex(infile, tabix)
to_outdir(ofile, outfile)
to_outdir(ofile + '.tbi', outidx)
