from os import path, makedirs
from glob import glob

{{runcmd}}
{{params2CmdArgs}}

# generate random regions
btparams = {}

btparams['l']    = {{args.len}}
btparams['n']    = 2*{{args.nperm}} # some regions do not have data
btparams['seed'] = {{i.seed}}
btparams['g']    = {{args.gsize | quote}}

randregfileUnsorted = path.join({{job.outdir | quote}}, 'randregions.bed.unsorted')
randregfile         = path.join({{job.outdir | quote}}, 'randregions.bed')
randregdir          = path.join({{job.outdir | quote}}, 'randregions')
cmd = '{{args.bedtools}} random %s > "%s"' % (params2CmdArgs(btparams, dash='-', equal=' '), randregfileUnsorted)
runcmd(cmd)

cmd = 'sort -k1,1 -k2,2n "%s" > "%s"' % (randregfileUnsorted, randregfile)
runcmd(cmd)

if not path.isdir(randregdir): makedirs(randregdir)
chroms = []
with open(randregfile) as f:
	prevchr = ''
	handler = None
	for line in f:
		if not line.strip(): continue
		chrom = line.split('\t')[0]
		if chrom != prevchr:
			chroms.append(chrom)
			if handler: handler.close()
			handler = open(path.join(randregdir, chrom + '.bed'), 'w')
			handler.write(line)
			prevchr = chrom
		else:
			handler.write(line)
	if handler:
		handler.close()


bwparams = {}
bwparams['skip-median'] = True

outfileUnsorted   = "{{o.outfile}}.unsorted"
outfileUnfiltered = "{{o.outfile}}.unfiltered"
for i, chrom in enumerate(chroms):
	regfile = path.join(randregdir, chrom + '.bed')
	bwfile  = glob(path.join({{args.consvdir | quote}}, chrom + '.*'))[0]
	cmd = '{{args.bwtool}} summary "%s" "%s" /dev/stdout %s | cat >> "%s"' % (regfile, bwfile, params2CmdArgs(bwparams, dash='-', equal=' '), outfileUnsorted)
	runcmd(cmd)

cmd = 'sort -k8nr "%s" > "%s"' % (outfileUnsorted, outfileUnfiltered)
runcmd(cmd)

# remove regions without data
with open(outfileUnfiltered) as f, open({{out.outfile | quote}}, 'w') as fout:
	for line in f:
		if '\tNA\tNA\tNA' in line: continue
		fout.write(line)

