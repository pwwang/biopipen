from os import path
from pyppl import Box
from bioprocs.utils.shell import cmdargs, runcmd
from bioprocs.utils.reference import vcfIndex

infile = {{i.infile | quote}}
outdir = {{o.outdir | quote}}
plink  = {{args.plink | quote}}
tabix  = {{args.tabix | quote}}
params = {{args.params | repr}}
infile = vcfIndex(infile, tabix = tabix)

params.vcf = infile
params['make-bed'] = True
params.out = path.join(outdir, {{i.infile | fn2 | quote}})

args = cmdargs(params, equal = ' ')
cmd = '{} {} 1>&2'.format(plink, args)
runcmd(cmd)