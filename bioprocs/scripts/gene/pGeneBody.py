import sys


{{ genenorm }}
{{ writeBedx }}

# get the genes
genes, _ = genenorm(
	{{in.infile | quote}}, 
	notfound = {{args.notfound | quote}},
	frm      = {{args.frm | quote}},
	to       = "genomic_pos_{{args.genome}},symbol",
	genome   = {{args.genome | quote}},
	tmpdir   = {{args.tmpdir | quote}},
	inopts   = {{args.inopts}},
	inmeta   = {{args.inmeta | lambda x: x if isinstance(x, list) or isinstance(x, dict) else '"' + x + '"'}}
)

writer = writeBedx({{out.outfile | quote}})
writer.meta.add(QUERY = 'The query gene name')
writer.writeHead()
for gene, hit in genes.items():
	pos      = hit['genomic_pos_{{args.genome}}']
	r        = readRecord()
	r.CHR    = 'chr' + str(pos['chr'])
	r.START  = pos['start']
	r.END    = pos['end']
	r.NAME   = hit['symbol']
	r.STRAND = '+' if pos['strand'] == 1 else '-'
	r.SCORE  = '0'
	r.QUERY  = gene
	writer.write(r)


