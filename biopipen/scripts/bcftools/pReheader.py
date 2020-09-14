import sys
from diot import Diot
from bioprocs.utils import shell2 as shell

infile   = {{i.infile | quote}}
hfile    = {{i.hfile | quote}}
samfile  = {{i.samfile | ?!:args.params.get('s', args.params.get('samples')) | repr}}
outfile  = {{o.outfile | quote}}
bcftools = {{args.bcftools | quote}}
params   = {{args.params | repr}}
nthread  = {{args.nthread | repr}}

shell.load_config(bcftools = bcftools)

params._    = infile
params.o    = outfile
if hfile:
	params.h = hfile
if samfile:
	if path.isfile(samfile):
		params.S = samfile
	else:
		params.s = samfile
cmd = shell.bcftools.reheader(**params).h.fg
sys.stderr.write("\n%s RUNNING %s\n%s\n%s\n\n" % (
	"-" * 40, "-" * 40, cmd.strcmd, "-" *  89))
cmd.run(wait=True)
