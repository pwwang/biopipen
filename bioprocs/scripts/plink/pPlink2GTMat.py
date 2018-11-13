
from os import path
from glob import glob
from bioprocs.utils import runcmd, cmdargs
from bioprocs.utils.tsvio import TsvReader, TsvWriter

# plink -bfile x --recode A-transpose --out x.txt 
# x.txt.traw

indir      = {{i.indir | quote}}
outfile    = {{o.outfile | quote}}
plink      = {{args.plink | quote}}
samid      = {{args.samid | quote}}
nors       = {{args.nors | quote}}
chroms     = {{args.chroms | repr}}

bedfile = glob(path.join(indir, '*.bed'))[0]
input   = path.splitext(bedfile)[0]
output  = path.splitext(outfile)[0]

params = {
	'bfile' : input,
	'recode': 'A-transpose',
	'out'   : output
}

cmd = '%s %s 1>&2' % (plink, cmdargs(params, equal = ' '))
runcmd(cmd)

fams = TsvReader(input + '.fam', ftype = 'nometa', delimit = ' ', head = False)
if samid == 'fid':
	header = "\t" + "\t".join(fams.dump(0)) + "\n"
elif samid == 'iid':
	header = "\t" + "\t".join(fams.dump(1)) + "\n"
else:
	header = "\t" + "\t".join(r[0] + '_' + r[1] for r in fams) + "\n"
fams.close()

gts = TsvReader(output + '.traw', ftype = 'nometa', skip = 1, head = False)

with open(outfile, 'w') as fout:
	fout.write(header)
	for gtline in gts:
		if gtline[0] in chroms:
			gtline[0] = chroms[gtline[0]]
		elif 'chr' + gtline[0] in chroms:
			gtline[0] = chroms[gtline[0]][3:] if chroms[gtline[0]].startswith('chr') else chroms[gtline[0]]
		if '_' in gtline[1]:
			gtline[1] = nors
		gtline[0] = '{chr}_{pos}_{rs}_{ref}_{alt}'.format(
			chr = gtline[0],
			pos = gtline[3],
			rs  = gtline[1],
			ref = gtline[4],
			alt = gtline[5]
		)
		del gtline[1:6]
		fout.write('\t'.join(gtline) + '\n')


