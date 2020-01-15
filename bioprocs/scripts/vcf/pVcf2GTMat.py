from collections import Counter
from pysam import VariantFile as Vcf
from bioprocs.utils.tsvio2 import TsvJoin

infile   = {{i.infile     | quote}}
outfile  = {{o.outfile    | quote}}
useid    = {{args.useid   | repr}}
na       = {{args.na      | quote}}
mingt    = {{args.mingt   | repr}}
novel    = {{args.novel   | quote}}
dbsnp    = {{args.dbsnp   | repr}}
bialt    = {{args.bialt   | repr}}
chrorder = {{args.chrorder.split(',') | repr}}
samname  = {{args.samname }}

chrorder = dict(zip(chrorder, list(range(len(chrorder)))))
if dbsnp:
	tsvjoin = TsvJoin(dbsnp, infile, cnames = False)
vcf = Vcf(infile)
samples = list(vcf.header.samples)
if callable(samname):
	samples = [samname(sam) for sam in samples]
sc = Counter(samples)
for s, c in sc.items():
	if c > 1:
		raise ValueError('Duplicated sample: {!r}'.format(s))

with open(outfile, 'w') as fout:
	fout.write("\t".join(samples) + "\n")
if mingt <= 1:
	mingt = len(samples) * mingt

def formatGT(gt):
	gt0 = gt[0]
	gt1 = gt[-1]
	try:
		gt0 = int(gt0)
		gt1 = int(gt1)
	except (ValueError, TypeError):
		return na
	gt = gt0 + gt1
	if bialt and gt > 2:
		return na
	return str(gt)

def formatRow(row, rs = novel):
	if bialt and row[3] not in list('AGTCagtc'):
		return None
	row[4] = row[4].rstrip(' ,')
	if bialt:
		row[4] = row[4].split(',')[0]
		if row[4] not in list('AGTCagtc'):
			return None
	gts = [formatGT(gt) for gt in row[9:]]
	if mingt > 0 and sum(1 for gt in gts if gt!=na) < mingt:
		return None
	if not useid or not row[2] or row[2] == '.':
		row[2] = rs
	if not rs: return None
	return "_".join(row[:5]) + "\t" + "\t".join(gts) + "\n"

if not dbsnp:
	reader = open(infile)
	writer = open(outfile, 'a')
	for line in reader:
		if line.startswith('#'):
			continue
		row = line.rstrip('\n').split('\t')
		row = formatRow(row)
		if row: writer.write(row)
	writer.close()
	reader.close()
else:
	match = lambda r1, r2: TsvJoin.compare(
		(chrorder[r1[0]], int(r1[1])),
		(chrorder[r2[0]], int(r2[1]))
	)
	def do(out, r1, r2):
		row = formatRow(r2, r1[2])
		if row: out.write(row)
	tsvjoin.join(do = do, match = match, outfile = outfile, outopts = dict(append = True))
