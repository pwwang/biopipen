import sys

{{ genenorm }}
# get the genes
genes, _ = genenorm(
	{{in.infile | quote}}, 
	col      = {{args.col}},
	notfound = {{args.notfound | quote}},
	frm      = {{args.frm | quote}},
	to       = "genomic_pos_{{args.genome}},symbol",
	header   = {{args.header}},
	genome   = {{args.genome | quote}},
	skip     = {{args.skip}},
	delimit  = {{args.delimit | quote}},
	tmpdir   = {{args.tmpdir | quote}},
	comment  = {{args.comment | quote}}
)

with open ({{out.outfile | quote}}, "w") as f:

	for gene, hit in genes.items():
		pos = hit['genomic_pos_{{args.genome}}']

		try:
			chr    = "chr" + str(pos['chr'])
			strand = pos['strand']
			start  = pos['start'] 
			end    = pos['end']
			f.write ("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chr, start, end, hit['symbol'], 0, ("+" if strand == 1 else "-"), gene))
		except TypeError:
			sys.stderr.write('Encounter TypeError, hit is: %s\n' % str(hit))
			raise

