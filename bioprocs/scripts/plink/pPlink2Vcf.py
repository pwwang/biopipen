from os import path, remove
from glob import glob
from bioprocs.utils import runcmd, cmdargs

indir   = {{i.indir | quote}}
outfile = {{o.outfile | quote}}
plink   = {{args.plink | quote}}
gz      = {{args.gz | repr}}
samid   = {{args.samid | quote}}
chroms  = {{args.chroms | repr}}

recode = ['vcf-fid' if samid == 'fid' else 'vcf-iid' if samid == 'iid' else 'vcf']
recode.append('tab')
if gz:
	outfile = outfile[:-3]

bedfile = glob(path.join(indir, '*.bed'))[0]
input   = path.splitext(bedfile)[0]
output  = path.splitext(outfile)[0]
outtmp  = output + '-tmp'

params = {
	'bfile' : input,
	'recode': recode,
	'out'   : outtmp
}
cmd = '%s %s 1>&2' % (plink, cmdargs(params, equal = ' '))
runcmd(cmd)

with open(outtmp + '.vcf', 'r') as fin, open(outfile, 'w') as fout:
	for line in fin:
		if line.startswith('#'):
			fout.write(line)
		else:
			items = line.split()
			if items[0] in chroms:
				items[0] = chroms[items[0]]
			elif 'chr' + items[0] in chroms:
				items[0] = chroms[items[0]][3:] if chroms[items[0]].startswith('chr') else chroms[items[0]]
			fout.write('\t'.join(items) + '\n')

remove(outtmp + '.vcf')
if gz:
	runcmd(['bgzip', outfile])

