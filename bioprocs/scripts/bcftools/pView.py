import sys
from os import path
from pyppl import Box
from bioprocs.utils import shell2 as shell
from bioprocs.utils.reference import vcfIndex

infile   = {{i.infile | quote}}
samfile  = {{i.samfile | quote}}
outfile  = {{o.outfile | quote}}
bcftools = {{args.bcftools | quote}}
params   = {{args.params | repr}}
gz       = {{args.gz | repr}}
nthread  = {{args.nthread | repr}}
tabix    = {{args.tabix | quote}}

vcfIndex(infile, tabix = tabix)

shell.load_config(bcftools = bcftools)

params._    = infile
params.o    = outfile
params.O    = 'z' if gz else 'v'
if samfile:
	if path.isfile(samfile):
		params.S = samfile
	else:
		params.s = samfile
cmd = shell.bcftools.view(**params).cmd
sys.stderr.write("\n%s RUNNING %s\n%s\n%s\n\n" % ("-" * 40, "-" * 40, cmd, "-" *  89))
