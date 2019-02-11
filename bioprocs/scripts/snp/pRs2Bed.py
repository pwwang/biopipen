from os import path
from pyppl import Box
from bioprocs.utils import shell
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter

snpfile   = {{ i.snpfile | quote}}
outfile   = {{ o.outfile | quote}}
inopts    = {{ args.inopts | repr}}
snpcol    = {{ args.snpcol or 0 | repr}}
dbsnp     = {{ args.dbsnp | quote}}
vcftools  = {{ args.vcftools | quote}}
sortby    = {{ args.sortby | quote}}
jobindir  = {{ job.indir | quote}}
joboutdir = {{ job.outdir | quote}}

shell.TOOLS.vcftools = vcftools
vcftools = shell.Shell(equal = ' ', dash = '--').vcftools

# snps
snplist = path.join(jobindir, path.basename(snpfile) + '.list')
reader  = TsvReader(snpfile, cnames = False)
writer  = TsvWriter(snplist)
for r in reader:
	writer.write([r[snpcol]])
reader.close()
writer.close()

params = Box()
params.snps = snplist
params.recode = True
params.out = path.join(joboutdir, 'tmp')
if dbsnp.endswith('.gz'):
	params.gzvcf = dbsnp
else:
	params.vcf = dbsnp
vcftools(**params).run()

reader = TsvReader(params.out + '.recode.vcf', cnames = False)
outfiletmp = outfile + '.tmp'
writer = TsvWriter(outfiletmp)
for r in reader:
	writer.write([r[0], r[1], r[1], r[2], 0, '+', r[3], r[4]])
reader.close()
writer.close()

if sortby == 'coord':
	shell.sort(k = ['1,1', '2,2n'], _ = outfiletmp, _stdout = outfile)
else:
	shell.sort(k = '4', _ = outfiletmp, _stdout = outfile)

shell.rm_rf(outfiletmp, params.out + '.recode.vcf')
