from pysam import VariantFile as Vcf
from bioprocs.utils.tsvio import TsvWriter, TsvRecord

infile  = {{i.infile   | quote}}
outfile = {{o.outfile | quote}}
rnames  = {{args.rnames | quote}}
na      = {{args.na     | quote}}
writer  = TsvWriter(outfile)

reader = Vcf(infile)
writer.meta.add('Variant', list(reader.header.samples))
writer.writeHead()
for r in reader.fetch():
	alts = r.alts
	if not alts: continue
	alts = list(alts)
	refalts = [r.ref] + alts
	for alt in alts:
		record = TsvRecord()
		if rnames == 'coord' or r.id == '.':
			record.Variant = '_'.join([r.chrom, str(r.pos), r.ref, alt])
		else:
			record.Variant = '_'.join([r.id, r.ref, alt])
		
		for samname, sample in r.samples.items():
			genotype = sample.get('GT', None)
			if not genotype or genotype == '.': 
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