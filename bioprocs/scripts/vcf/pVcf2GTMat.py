import tabix
from pysam import VariantFile as Vcf
from bioprocs.utils.tsvio import TsvWriter, TsvRecord

infile  = {{i.infile     | quote}}
outfile = {{o.outfile    | quote}}
useid   = {{args.useid   | R}}
na      = {{args.na      | quote}}
novel   = {{args.novel   | quote}}
dbsnp   = {{args.dbsnp   | quote}}
bialt   = {{args.bialt   | R}}
writer  = TsvWriter(outfile)

tb = tabix.open(dbsnp) if dbsnp else None
def getRsName(chr, pos):
	if not tb:
		return novel
	try:
		records = list(tb.query(chr, pos-1, pos))
		if not records:
			return novel
		return records[0][2]
	except tabix.TabixError:
		return novel

reader = Vcf(infile)
writer.meta.add('Variant')
writer.meta.add(*tuple(reader.header.samples))
writer.writeHead()
for r in reader.fetch():
	alts = r.alts
	if not alts: continue
	if bialt and len(alts) != 1: continue
	if bialt and len(r.ref) != 1: continue
	alts    = list(alts)
	refalts = [r.ref] + alts
	name    = r.id if useid and r.id and r.id != '.' else getRsName(r.chrom, r.pos)
	for alt in alts:
		record = TsvRecord()
		record.Variant = '{chr}_{pos}_{name}_{ref}_{alt}'.format(
			chr  = r.chrom,
			pos  = r.pos,
			name = name,
			ref  = r.ref,
			alt  = alt
		)
		
		for samname, sample in r.samples.items():
			genotype = sample.get('GT', None)
			if not genotype or None in genotype or '.' in genotype: 
				genotype = na
			else:
				gtidx0 = int(genotype[0])
				gtidx1 = int(genotype[-1])
				gt0    = refalts[gtidx0]
				gt1    = refalts[gtidx1]
				if gt0+gt1 == r.ref+r.ref:
					genotype = 0
				elif gt0+gt1 == r.ref+alt or gt0+gt1 == alt+r.ref:
					genotype = 1
				elif gt0+gt1 == alt+alt:
					genotype = 2
				else:
					genotype = na
					
			record[samname] = genotype
		writer.write(record)
			
writer.close()