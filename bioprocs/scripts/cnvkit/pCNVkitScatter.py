from os import path
from copy import deepcopy
from pyppl import Box
from bioprocs.utils import runcmd, cmdargs

cnvkit   = {{args.cnvkit | quote}}
cnrfile  = {{i.cnrfile | quote}}
cnsfile  = {{i.cnsfile | quote}}
outdir   = {{o.outdir | quote}}
regions  = {{args.regions | repr}}
nthread  = {{args.nthread | repr}}
params   = {{args.params}}

openblas_nthr = "export OPENBLAS_NUM_THREADS={nthread}; export OMP_NUM_THREADS={nthread}; export NUMEXPR_NUM_THREADS={nthread}; export MKL_NUM_THREADS={nthread}".format(nthread = nthread)

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

	cmd = '{openblas}; {cnvkit} scatter \'{cnrfile}\' {params}'
	runcmd(cmd.format(**Box(
		openblas = openblas_nthr,
		cnvkit   = cnvkit,
		params   = cmdargs(iparams, equal = ' '),
		cnrfile  = cnrfile
	)))


