from os import path
from copy import deepcopy
from pyppl import Box
from bioprocs.utils import shell

cnvkit   = {{args.cnvkit | quote}}
cnrfile  = {{i.cnrfile | quote}}
cnsfile  = {{i.cnsfile | quote}}
outdir   = {{o.outdir | quote}}
regions  = {{args.regions | repr}}
nthread  = {{args.nthread | repr}}
params   = {{args.params}}

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
		genes = ''
	else:
		parts = region.split(':')
		if len(parts) == 1:
			rgion = parts[0]
		elif len(parts) == 2:
			import re
			if re.match(r'^\d+\-\d+$', parts[1]):
				rgion = ':'.join(parts)
				genes = ''
			else:
				rgion = parts[0]
				genes = parts[1]
		else:
			rgion = ':'.join(parts[:2])
			genes = parts[2]

	stem = region.replace(':', '_').replace(',', '_')
	if not stem: stem = 'whole-genome'

	outfile   = path.join(outdir, stem + '.pdf')
	iparams   = deepcopy(params)
	iparams.o = outfile
	
	if cnsfile:
		iparams.s = cnsfile

	if rgion:
		iparams.c = rgion
	
	if genes:
		iparams.g = genes

	ckshell.scatter(cnrfile, **iparams).run()


