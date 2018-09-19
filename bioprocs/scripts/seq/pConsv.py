from os import path, makedirs
from glob import glob

{{runcmd}}

# split input bed file
beddir = path.join({{job.outdir | quote}}, 'chrom-regions')
if not path.isdir(beddir): makedirs(beddir)

bedfile = path.join({{job.outdir | quote}}, "{{i.bedfile | bn}}.sorted")
cmd = 'grep -v "^#" {{i.bedfile | quote}} | cut -f1,2,3,4,5,6  | sort -k1,1 -k2,2n > "%s"' % bedfile
runcmd(cmd)

chroms = []
with open(bedfile) as f:
	prevchr = ''
	handler = None
	for line in f:
		if not line.strip(): continue
		chrom = line.split('\t')[0]
		if chrom != prevchr:
			chroms.append(chrom)
			if handler: handler.close()
			handler = open(path.join(beddir, chrom + '.bed'), 'w')
			handler.write(line)
			prevchr = chrom
		else:
			handler.write(line)
	if handler:
		handler.close()

{% if args.pval %}
def getPval(num, dist):
	for i, d in enumerate(dist):
		if d < num: 
			return float(i)/float(len(dist))
	return 1.0

with open({{i.permfile | quote}}) as f:
	nulldist = [float(line.split('\t')[7].strip()) for line in f.read().splitlines() if line.strip()]
{% endif %}

# handle each chrom
outfileNoPval = "{{o.outfile}}.nopval"
for chrom in chroms:
	regfile = path.join(beddir, chrom + '.bed')
	bwfile  = glob(path.join({{args.consvdir | quote}}, chrom + '.*'))[0]
	cmd = '{{args.bwtool}} summary "%s" "%s" /dev/stdout -skip-median | cat >> "%s"' % (regfile, bwfile, outfileNoPval)
	runcmd (cmd)

outs = {}
with open(outfileNoPval) as f:
	for line in f:
		line  = line.strip()
		if not line: continue
		parts = line.split('\t')
		outs['\t'.join(parts[:3])] = 0 if parts[7] == 'NA' else float(parts[7])

with open({{i.bedfile | quote}}) as f, open({{out.outfile | quote}}, 'w') as fout:
	for line in f:
		line  = line.strip()
		if not line: continue
		if line.startswith('#'):
			{% if args.pval %}
			fout.write(line + '\tConservation\tConsvPval\n')
			{% else %}
			fout.write(line + '\tConservation\n')
			{% endif %}
			continue
		parts = line.split('\t')
		consv = outs['\t'.join(parts[:3])]
		{% if args.pval %}
		pval  = getPval(consv, nulldist)
		if pval >= {{args.pval}}: continue
		fout.write(line + '\t%.3f\t%.2E\n' % (consv, pval))
		{% else %}
		fout.write(line + '\t%.3f\n' % (consv))
		{% endif %}

