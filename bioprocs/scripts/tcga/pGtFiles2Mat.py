from os import path
from bioprocs.utils import logger

infiles = {{in.infiles | repr}}
outfile = {{out.outfile | quote}}
rsmap   = {{args.rsmap | : None if not v1 else v1 | repr}}
fn2sam  = {{args.fn2sam}}

def readone(infile, getsnps = True):
	logger.info('Reading ' + infile + ' ...')
	gts        = []
	snps       = []
	gts_append = gts.append
	if getsnps:
		snps_append = snps.append
	with open(infile) as fin:
		fin.readline()
		fin.readline()
		for line in fin:
			line = line.strip()
			if not line: continue
			snp, gt = line.split()[:2]
			if getsnps:
				snps_append(snp)
			gts_append(gt)
	return snps, gts

ret = None
samples = ['Variant']
for i, infile in enumerate(infiles):
	sam = path.splitext(path.basename(infile))[0]
	if fn2sam:
		samples.append(fn2sam(sam))
	else:
		samples.append(sam)
	snps, gts = readone(infile, i == 0)
	if snps: 
		ret = [(snps[i], [gt]) for i, gt in enumerate(gts)]
	else:
		[ret[i][1].append(gt) for i, gt in enumerate(gts)]
		
if rsmap:
	logger.info('Reading rs mapping file ...')
	with open(rsmap) as fin:
		snp2rs = {snp:rs for line in fin if line.strip() for rs, snp in [line.split()[:2]]}
	ret = [(snp2rs[r[0]], r[1]) for r in ret if snp2rs[r[0]].startswith('rs')]

logger.info('Saving results ...')
with open(outfile, 'w') as fout:
	fout.write('\t'.join(samples) + '\n')
	for rs, gts in ret:
		fout.write(rs + '\t' + '\t'.join(gts) + '\n')
