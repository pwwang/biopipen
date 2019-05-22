from os import path
from pyppl import Box
from bioprocs.utils import ensureBox, shell2 as shell
from bioprocs.utils.reference import vcfIndex

infile = {{i.infile | quote}}
outdir = {{o.outdir | quote}}
plink  = {{args.plink | quote}}
tabix  = {{args.tabix | quote}}
params = {{args.params}}
params = ensureBox(params)

# resolve plink 1.x --set-missing-var-ids doesn't distinguish $1, $2,... for ref and alts
if 'set-missing-var-ids' in params and "$" in params['set-missing-var-ids']:
	tmpfile = path.join(outdir, 'withvarids.vcf')
	import vcf
	reader = vcf.Reader(filename = infile)
	stream = open(tmpfile, 'w')
	writer = vcf.Writer(stream, reader)
	for r in reader:
		if not r.ID:
			all_alts = [r.REF] + list(r.ALT)
			r.ID     = params['set-missing-var-ids'].replace(
				'@', r.CHROM).replace('#', str(r.POS))
			for i, a in enumerate(all_alts):
				r.ID = r.ID.replace('${}'.format(i+1), str(a))
		writer.write_record(r)
	stream.close()
	infile = tmpfile
infile = vcfIndex(infile, tabix = tabix)

params.vcf = infile
params['make-bed'] = True
params.out = path.join(outdir, {{i.infile | fn2 | quote}})

shell.fg.plink(**params)
