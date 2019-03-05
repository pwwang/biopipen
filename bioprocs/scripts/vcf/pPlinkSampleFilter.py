from os import path
from glob import glob
from pyppl import Box
from bioprocs.utils import shell

indir   = {{i.indir | quote}}
samfile = {{i.samfile | quote}}
outdir  = {{o.outdir | quote}}
plink   = {{args.plink | quote}}
keep    = {{args.keep | repr}}
fam     = {{args.fam | repr}}
params  = {{args.params | repr}}

plink = shell.Shell(dict(plink = plink), equal = ' ').plink

bedfile = glob(path.join(indir, '*.bed'))[0]
input   = bedfile[:-4]
output  = path.join(outdir, path.basename(input))

params.bfile   = bedfile
params.out     = output
params.nthread = nthread
params['make-bed'] = True

key = 'keep' if keep and not fam else 'keep-fam' if keep and fam else 'remove' if not keep and not fam else 'remove-fam'
params[key] = samfile

plink(**params).run()
