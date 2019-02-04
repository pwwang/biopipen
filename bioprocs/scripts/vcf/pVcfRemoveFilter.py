from pysam import VariantFile as Vcf
from pyppl import Box
from bioprocs.utils import alwaysList

infile   = {{ i.infile | quote}}
outfile  = {{ o.outfile | quote}}
rmfilter = {{ args.rmfilter | repr}}
if rmfilter:
	rmfilter = alwaysList(rmfilter)

invcf  = Vcf(infile)
outvcf = open(outfile, 'w')
outvcf.write(str(invcf.header))
for rec in invcf.fetch():
	parts   = str(rec).split('\t')
	filters = parts[6].split(';')
	if not rmfilter:
		filters = 'PASS'
	else:
		filters = ';'.join(f for f in filters if f not in rmfilter)
		filters = filters or 'PASS'
	parts[6] = filters
	outvcf.write('\t'.join(parts))
outvcf.close()


