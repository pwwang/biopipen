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
genes = list(sorted(set(reader.dump('Hugo_Symbol'))))
reader.rewind()
samples = list(sorted(set(reader.dump('Tumor_Sample_Barcode'))))
reader.rewind()

matrix = {gene: {sample: 0 for sample in samples} for gene in genes}
for r in reader:
	if r.Variant_Classification not in classes:
		continue
	if binary:
		matrix[r.Hugo_Symbol][r.Tumor_Sample_Barcode] = 1
	else:
		matrix[r.Hugo_Symbol][r.Tumor_Sample_Barcode] += 1

writer = TsvWriter(outfile)
writer.cnames = ['Gene'] + samples
writer.writeHead()

for gene in genes:
	counts = [matrix[gene][sample] for sample in samples]
	if sum(counts) == 0 and not inc_all:
		continue
	writer.write([gene] + counts)
writer.close()
