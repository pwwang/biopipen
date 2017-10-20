import sys
from mygene import MyGeneInfo

genome2species = {
	'hg19': 'human',
	'hg38': 'human',
	'mm9' : 'mouse',
	'mm10': 'mouse',
}
species = genome2species[{{args.genome | quote}}]

mg = MyGeneInfo()
genes = [line.split()[0] for line in open({{in.genefile | quote}}) if line.strip()]
genes = mg.querymany(genes, size=1, fields="genomic_pos_{{args.genome}},symbol", scopes="symbol,alias", species=species)

with open ("{{out.outfile}}", "w") as f:
	for hit in genes:
		if not 'genomic_pos_{{args.genome}}' in hit:
			sys.stderr.write('pyppl.log.warning: Cannot find position for gene: %s\n' % hit['query'])
			continue
		
		pos = hit['genomic_pos_{{args.genome}}']
		
		try:
			chr    = "chr" + str(pos['chr'])
			strand = pos['strand']
			tss    = pos['start'] if strand == 1 else pos['end']
			pstart = tss - {{args.up}}
			pend   = tss + {{args.down}}
			f.write ("%s\t%s\t%s\t%s\t%s\t%s\n" % (chr, pstart, pend, hit['symbol'], 0, ("+" if strand == 1 else "-")))
		except TypeError:
			sys.stderr.write('Encounter TypeError, hit is: %s\n' % str(hit))