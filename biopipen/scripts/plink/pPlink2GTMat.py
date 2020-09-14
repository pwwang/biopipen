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
outsnp = {{o.outsnp | quote}}
plink = {{args.plink | quote}}
samid = {{args.samid | quote}}
snpid = {{args.snpid | quote}}
addchr = {{args.addchr | repr}}
nors = {{args.nors | quote}}
refallele = {{args.refallele | repr}}
chroms = {{args.chroms | repr}}
bcftools = {{args.bcftools | quote}}

bedfile = glob(path.join(indir, '*.bed'))[0]
prefix = path.splitext(bedfile)[0]
output = path.splitext(outfile)[0]

params = {
    'bfile': prefix,
    'recode': 'A-transpose',
    'out': output
}

shell.load_config(plink=plink, bcftools=bcftools)

shell.plink(**params).fg

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
        snpname = gtline[1]
    else:
        if gtline[0] in chroms:
            gtline[0] = chroms[gtline[0]]
        elif 'chr' + gtline[0] in chroms:
            gtline[0] = chroms['chr' + gtline[0]]
        if addchr and not gtline[0].startswith('chr'):
            gtline[0] = 'chr' + gtline[0]
        snpname = snpid.format(
            chr=gtline[0],
            pos=gtline[3],
            rs=gtline[1] if 'rs' in gtline[1] else nors,
            ref=gtline[4], # NOT reliable!!
            alt=gtline[5] # NOT reliable!!
        )
    gtline[0] = snpname
    del gtline[1:6]
    writer.write(list(gtline.values()))
writer.close()

# plink 1.9 does not keep track of reference allele
# we need either a vcf with reference alleles or
# a bed6+ file with ref/alt alleles as 7th and 8th columns
if refallele:
    params = {
        'bfile': prefix,
        'recode': 'vcf-iid',
        'out': outsnp,
    }
    if refallele.endswith('.vcf.gz') or refallele.endswith('.vcf'):
        # get the vcf file first to extract information
        params['output-chr'] = 'chr' + chroms["26"]
        shell.plink(**params).fg
        shell.bcftools.query(
            R=outsnp + '.vcf',
            o=outsnp,
            f='%CHROM\\t%POS0\\t%END\\t%ID\\t0\\t+\\t%REF\\t%ALT{0}\\n',
            _=refallele
        ).fg
        refallele = outsnp

    params['a2-allele'] = [refallele, 7, 4, '#']
    shell.plink(**params).fg

    snpreader = TsvReader(outsnp + '.vcf', cnames=False)
    snpwriter = TsvWriter(outsnp)
    for snp in snpreader:
        snpwriter.write([
            ('chr' + snp[0]
             if addchr and not snp[0].startswith('chr')
             else snp[0]),
            snp[1],
            snp[1],
            snp[2],
            0,
            '+', snp[3], snp[4]
        ])
    snpwriter.close()
    snpreader.close()
else:
    open(outsnp, 'w').close()
