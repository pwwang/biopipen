from bioprocs.utils.tsvio2 import TsvReader, TsvWriter
from bioprocs.utils.constants import tcgamaf

{% from pyppl.utils import always_list %}
infile  = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
classes = {{args.classes | always_list | repr}}
binary  = {{args.binary | bool}}
inc_all = {{args.all | bool}}

for group, clses in tcgamaf.classes.items():
	if group in classes:
		classes.remove(group)
		classes.extend(clses)

reader = TsvReader(infile)
mutations = [
	'%s:%s' % (chrom, start)
	for chrom, start in sorted(set(reader.dump(['Chromosome', 'Start_Position'])))
]
reader.rewind()
samples = list(sorted(set(reader.dump('Tumor_Sample_Barcode'))))
reader.rewind()

matrix = {mutation: {sample: 0 for sample in samples} for mutation in mutations}
for r in reader:
	if r.Variant_Classification not in classes:
		continue
	key = '%s:%s' % (r.Chromosome, r.Start_Position)
	if binary:
		matrix[key][r.Tumor_Sample_Barcode] = 1
	else:
		matrix[key][r.Tumor_Sample_Barcode] += 1

writer = TsvWriter(outfile)
writer.cnames = ['Mutation'] + samples
writer.writeHead()

for mutation in mutations:
	counts = [matrix[mutation][sample] for sample in samples]
	if sum(counts) == 0 and not inc_all:
		continue
	writer.write([mutation] + counts)
writer.close()
