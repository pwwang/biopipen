"""Script for plink.pPlink2GTMat"""
# pylint: disable=undefined-variable,unused-import,invalid-name
# pylint: disable=unsupported-assignment-operation,not-a-mapping
# pylint: disable=not-an-iterable,unsubscriptable-object

from os import path
from glob import glob
from diot import Diot
from bioprocs.utils import shell2 as shell
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter

# plink -bfile x --recode A-transpose --out x.txt
# x.txt.traw

indir = {{i.indir | quote}}
outfile = {{o.outfile | quote}}
plink = {{args.plink | quote}}
samid = {{args.samid | quote}}
snpid = {{args.snpid | quote}}
addchr = {{args.addchr | repr}}
nors = {{args.nors | quote}}
chroms = {{args.chroms | repr}}

bedfile = glob(path.join(indir, '*.bed'))[0]
prefix = path.splitext(bedfile)[0]
output = path.splitext(outfile)[0]

params = {
    'bfile': prefix,
    'recode': 'A-transpose',
    'out': output
}

shell.load_config(plink=plink)

shell.fg.plink(**params)

fams = TsvReader(prefix + '.fam', delimit=' ', cnames=False)
if samid == 'fid':
    header = ['ID'] + fams.dump(0)
elif samid == 'iid':
    header = ['ID'] + fams.dump(1)
else:
    header = ['ID'] + [r[0] + '_' + r[1] for r in fams]
fams.close()

gts = TsvReader(output + '.traw', skip=1, cnames=False)
writer = TsvWriter(outfile)
writer.cnames = header
writer.writeHead()

for gtline in gts:
    if snpid == 'raw':
        gtline[0] = gtline[1]
    else:
        if gtline[0] in chroms:
            gtline[0] = chroms[gtline[0]]
        elif 'chr' + gtline[0] in chroms:
            gtline[0] = chroms['chr' + gtline[0]]
        if addchr and not gtline[0].startswith('chr'):
            gtline[0] = 'chr' + gtline[0]
        gtline[0] = snpid.format(
            chr=gtline[0],
            pos=gtline[3],
            rs=gtline[1] if 'rs' in gtline[1] else nors,
            ref=gtline[4],
            alt=gtline[5]
        )
    del gtline[1:6]
    writer.write(list(gtline.values()))
writer.close()
