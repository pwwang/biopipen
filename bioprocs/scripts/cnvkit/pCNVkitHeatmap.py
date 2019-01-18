from os import path
from copy import deepcopy
from pyppl import Box
from bioprocs.utils import runcmd, cmdargs

cnvkit   = {{args.cnvkit | quote}}
cnfiles  = {{i.cnfiles | asquote | quote}}
outdir   = {{o.outdir | quote}}
regions  = {{args.regions | repr}}
params   = {{args.params}}
nthread  = {{args.nthread | repr}}

for region in regions:
	if not region:
		rgion = ''
	else:
		parts = region.split(':')
		if len(parts) == 1:
			rgion = parts[0]
		elif len(parts) == 2:
			import re
			if re.match(r'^\d+\-\d+$', parts[1]):
				rgion = ':'.join(parts)
			else:
				rgion = parts[0]
		else:
			rgion = ':'.join(parts[:2])

	stem = region.replace(':', '_').replace(',', '_')
	if not stem: stem = 'whole-genome'

	outfile   = path.join(outdir, stem + '.pdf')
	iparams   = deepcopy(params)
	iparams.o = outfile

	if rgion:
		iparams.c = rgion

	cmd = '{nthr}; {cnvkit} heatmap {cnfiles} {params}'
	runcmd(cmd.format(**Box(
		nthr     = "export OPENBLAS_NUM_THREADS={nthread}; export OMP_NUM_THREADS={nthread}; export NUMEXPR_NUM_THREADS={nthread}; export MKL_NUM_THREADS={nthread}".format(nthread = nthread),
		cnvkit   = cnvkit,
		params   = cmdargs(iparams, equal = ' '),
		cnfiles  = cnfiles
	)))


