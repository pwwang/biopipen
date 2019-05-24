from os import path
from pyppl import Box
from bioprocs.utils import shell2 as shell

cnvkit   = {{args.cnvkit | quote}}
cnfiles  = {{i.cnfiles | repr}}
outdir   = {{o.outdir | quote}}
regions  = {{args.regions | list}}
params   = {{args.params | repr}}
nthread  = {{args.nthread | repr}}

shell.load_config(cnvkit = dict(
	_exe = cnvkit,
	_env = dict(
		OPENBLAS_NUM_THREADS = str(nthread),
		OMP_NUM_THREADS      = str(nthread),
		NUMEXPR_NUM_THREADS  = str(nthread),
		MKL_NUM_THREADS      = str(nthread)
	),
	_cwd = outdir
))

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
	iparams   = params.copy()
	iparams.o = outfile

	if rgion:
		iparams.c = rgion

	shell.fg.cnvkit.heatmap(*cnfiles, **iparams)

