from os import path
from pyppl import Box
from vcf import Reader as Vcf
from bioprocs.utils import runcmd, cmdargs, logger

infile  = {{i.infile | quote}}
idxfile = infile + '.tbi'
if not path.isfile(idxfile):
	raise ValueError('Vcf file needs to be indexed')

outdir = {{o.outdir | quote}}
plink  = {{args.plink | repr}}
params = {{args.params | repr}}

params.vcf = infile
params['make-bed'] = True
params.out = path.join(outdir, {{i.infile | fn2 | quote}})

args = cmdargs(params, equal = ' ')
cmd = '{} {} 1>&2'.format(plink, args)
runcmd(cmd)