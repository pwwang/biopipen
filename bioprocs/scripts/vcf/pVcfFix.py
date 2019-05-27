import gzip
from pysam import VariantFile # used to operate header
from pyppl import Box
from bioprocs.utils.reference import vcfIndex

infile  = {{i.infile | quote}}
outfile = {{o.outfile | quote}}
fixes   = {{args.fixes | repr}}
tabix   = {{args.tabix | quote}}

vcfIndex(infile, tabix = tabix)

# opoerator header

openfun = gzip.open if infile.endswith('.gz') else open
with openfun(infile, 'rt', errors='replace') as fin, \
	openfun(outfile, 'at', errors='replace') as fout:
	for line in fin:
		if line.startswith('#'):
			# if no header operated
			fout.write(line)
			continue
		parts = line.strip().split('\t')
		if fixes.addChr:
			parts[0] = parts[0] if parts[0].startswith('chr') else 'chr' + parts[0]
		info = parts[7]
		if fixes.clinvarLink:
			info = ';'.join(inf for inf in info.split(';') if not inf.startswith('<a href'))
		parts[7] = info

		fout.write('\t'.join(parts) + '\n')

