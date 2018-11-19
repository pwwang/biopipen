import vcf
from bioprocs.utils.tsvio import TsvWriter, TsvRecord

infile  = {{i.infile | quote}}
outfile = {{o.outfile | quote}}

reader = vcf.Reader(filename = infile)
writer = TsvWriter(outfile)
writer.meta.add(
	'mutation_id',
	'ref_counts',
	'var_counts',
	'normal_cn',
	'minor_cn',
	'major_cn',
	'variant_case',
	'variant_freq',
	'genotype'
)
while True:
	try:
		r = next(reader)
		record = TsvRecord()
		record.normal_cn = 2
		# mutation_id	ref_counts	var_counts	normal_cn	minor_cn	major_cn	variant_case	variant_freq	genotype
		# NA12156:BB:chr2:175263063	3812	14	2	0	2	NA12156	0.0036591741	BB
		for samdata in record.samples:
			

	except StopIteration:
		break
	except:
		continue
