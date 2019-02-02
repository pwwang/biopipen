from os import path
from copy import deepcopy
from pyppl import Box
from bioprocs.utils import shell

cnvkit   = {{args.cnvkit | quote}}
cnfiles  = {{i.cnfiles | repr}}
outdir   = {{o.outdir | quote}}
regions  = {{args.regions | repr}}
params   = {{args.params}}
nthread  = {{args.nthread | repr}}

shell.TOOLS['cnvkit'] = cnvkit
envs = dict(
	OPENBLAS_NUM_THREADS = nthread,
	OMP_NUM_THREADS      = nthread,
	NUMEXPR_NUM_THREADS  = nthread,
	MKL_NUM_THREADS      = nthread
)
ckshell = shell.Shell(subcmd = True, equal = ' ', envs = envs, cwd = outdir).cnvkit

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

	ckshell.heatmap(*cnfiles, **iparams)

