from os import path
from glob import glob
from pyppl import Box
from bioprocs.utils import shell

indir   = {{i.indir | quote}}
samfile = {{i.samfile | quote}}
outdir  = {{o.outdir | quote}}
plink   = {{args.plink | quote}}
nthread = {{args.nthread | repr}}
keep    = {{args.keep | repr}}
fam     = {{args.fam | repr}}
samname = {{args.samid | repr}}
params  = {{args.params | repr}}

plink = shell.Shell(dict(plink = plink), equal = ' ').plink

bedfile = glob(path.join(indir, '*.bed'))[0]
input   = bedfile[:-4]
output  = path.join(outdir, path.basename(input))

# prepare sample file
if samname == 'both':
	pass
else:
	samplefile = path.join(outdir, 'samples.txt')
	with open(samfile, 'r') as fin, open(samplefile, 'w') as fout:
		for line in fin:
			fout.write("\t".join(line.strip().split("\t")[:1] * 2) + "\n")
	samfile = samplefile

params.bfile   = input
params.out     = output
params.nthread = nthread
params['make-bed'] = True

key = 'keep' if keep and fam else 'keep-fam' if keep and not fam else 'remove' if not keep and not fam else 'remove-fam'
params[key] = samfile

plink(**params).run()
