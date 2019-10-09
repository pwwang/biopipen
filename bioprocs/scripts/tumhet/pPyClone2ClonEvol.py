
# outfile should be:
"""
cluster	gene	is.driver	P.vaf	R.vaf	P.ccf	R.ccf	P.ref.count	P.var.count	P.depth	R.ref.count	R.var.count	R.depth	P	R
1	SMC3	TRUE	50.28	49.13	100.56	98.26	3027	2611	5639	2035	1990	4024	50.28	49.13
1	PTPRT	TRUE	45.24	45.31	90.48	90.62	2159	1799	3958	1715	1579	3294	45.24	45.31
1	-	FALSE	43.83	41.37	87.66	82.74	2931	2364	5296	2224	1501	3725	43.83	41.37
"""

# pyclone loci file
"""
mutation_id	sample_id	cluster_id	cellular_prevalence	cellular_prevalence_std	variant_allele_frequency
chr10:17875816	NA12878	2	0.828865827591	0.0355432431918	0.904883227176
chr10:17875816	NA12156	2	0.741264908504	0.11050033704	0.560311284047
chr10:17875816	NA19240	2	0.837320950405	0.0572230410102	0.89578313253
"""

# ref/alt counts from pyclone's yaml file:
"""
mutations:
- id: chr10:91143374
  ref_counts: 5506
  states:
  - g_n: AA
    g_r: AA
    g_v: AB
    prior_weight: 1
  - g_n: AA
    g_r: AA
    g_v: BB
    prior_weight: 1
  var_counts: 33
"""
from pathlib import Path
import yaml
from bioprocs.utils.tsvio2 import TsvReader, TsvWriter, TsvRecord
from bioprocs.utils import shell2 as shell
indir    = Path({{i.indir | quote}})
outfile  = Path({{o.outfile | quote}})
refgene  = {{args.refgene | quote}}
bedtools = {{args.bedtools | quote}}
drivers  = {{args.drivers | str}}
shell.load_config(bedtools = bedtools)

# get samples
samples = [path.stem for path in (indir/'yaml').glob('*.yaml')]

# get basic mutation information from loci file
locifile = indir / 'tables' / 'loci.tsv'
locireader = TsvReader(locifile)
records = {}
for loci in locireader:
	if loci.mutation_id not in records:
		records[loci.mutation_id] = TsvRecord()
	records[loci.mutation_id][loci.sample_id + '.vaf'] = loci.variant_allele_frequency
	records[loci.mutation_id][loci.sample_id] = loci.variant_allele_frequency
	# is it ok to use cellular_prevalence as ccf?
	records[loci.mutation_id][loci.sample_id + '.ccf'] = loci.cellular_prevalence
	records[loci.mutation_id]['cluster'] = int(loci.cluster_id) + 1

# get the counts from yaml file
for ymlfile in (indir/'yaml').glob('*.yaml'):
	with ymlfile.open() as ymlstrm:
		mutations = yaml.safe_load(ymlstrm)
		for mutation in mutations['mutations']:
			if mutation['id'] not in records:
				records[mutation['id']] = TsvRecord()
			records[mutation['id']][ymlfile.stem + '.ref.count'] = mutation['ref_counts']
			records[mutation['id']][ymlfile.stem + '.var.count'] = mutation['var_counts']
			records[mutation['id']][ymlfile.stem + '.depth'] = int(mutation['ref_counts']) + int(mutation['var_counts'])

# annotate mutations with genes
# construct mutation bed file
bedwriter = TsvWriter(outfile.parent / 'mutations.bed')
for mutation in records:
	chrom, start = mutation.split(':')
	end = int(start) + 1
	bedwriter.write([chrom, start, end, mutation, 0, '+'])
bedwriter.close()

mutations_with_genes = outfile.parent / 'mutations_with_genes.txt'
shell.bedtools.intersect(a = bedwriter.filename, b = refgene, wo = True, _out = str(mutations_with_genes))

bedreader = TsvReader(mutations_with_genes)
for bedr in bedreader:
	if bedr[3] not in records:
		continue
	# gene_id "IFIT1B"; transcript_id "IFIT1B"; gene_name "IFIT1B";
	gene = bedr[14].split('; ')[0][8:].strip('"')
	if 'gene' in records[bedr[3]]:
		records[bedr[3]].gene += ',' + gene
	else:
		records[bedr[3]].gene = gene
	if gene in drivers:
		records[bedr[3]]["is.driver"] = 'TRUE'

writer = TsvWriter(outfile)
writer.cnames = ['cluster', 'gene', 'is.driver']
for suffix in ('.vaf', '.ccf', '.ref.count', '.var.count', '.depth', ''):
	for sample in samples:
		writer.cnames.append(sample + suffix)
writer.writeHead()
# fill the defaults
for mutation, record in records.items():
	if 'cluster' not in record:
		#record['cluster'] = 'NA'
		# skip mutations with cluster
		continue
	if 'gene' not in record:
		record['gene'] = '-'
	if 'is.driver' not in record:
		record['is.driver'] = 'FALSE'
	for sample in samples:
		for suffix in ('.vaf', '.ccf', '.ref.count', '.var.count', '.depth', ''):
			if sample + suffix not in record:
				record[sample + suffix] = 'NA'
	writer.write(record)
writer.close()
