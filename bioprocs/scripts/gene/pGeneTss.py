import sys, json
from os import makedirs, path
from hashlib import md5
from mygene import MyGeneInfo

{% if args.header %}
skip = {{args.skip}} + 1
{% else %}
skip = {{args.skip}}
{% endif %}

# get the genes
delimit = {{args.delimit | quote}}
comment = {{args.comment | quote}}
tmpdir  = {{args.tmpdir | quote}}
col     = {{args.col}}
with open({{in.infile | quote}}) as f:
	genes = list(sorted(set([line.split(delimit)[col] for i, line in enumerate(f.read().splitlines()) if i >= skip and line.strip() and not line.startswith(comment)])))

genome2species = {
	'hg19': 'human',
	'hg38': 'human',
	'mm9' : 'mouse',
	'mm10': 'mouse',
}
species = genome2species[{{args.genome | quote}}]

frm = {{args.frm | quote}}
to  = "genomic_pos_{{args.genome}},symbol"
cachedir = path.join(tmpdir, 'mygene-cache')
if not path.isdir(cachedir): makedirs(cachedir)

uid   = md5(''.join(genes) + str(frm) + str(to) + str(species)).hexdigest()
cache = path.join(cachedir, 'mygeneinfo.%s' % uid)
if path.isfile(cache):
	with open(cache) as f: mgret = json.load(f)
else:
	mgret = MyGeneInfo().getgenes(genes, scopes = frm, fields = to, species = species)
	with open(cache, 'w') as f:	json.dump(mgret, f)

hits = {}
for hit in mgret:
	if not 'genomic_pos_{{args.genome}}' in hit:
		continue
	
	q = hit['query']
	if not q in hits: hits[q] = []
	hits[q].append(hit)

with open ("{{out.outfile}}", "w") as f:

	for gene in genes:
		if not gene in hits:
			sys.stderr.write('pyppl.log.warning: Cannot find coordinates for gene: %s\n' % gene)
			continue
		
		# multiple hits, select the one with greatest score
		hit = sorted(hits[gene], key=lambda x: x['_score'], reverse=True)[0]
		pos = hit['genomic_pos_{{args.genome}}']

		try:
			chr    = "chr" + str(pos['chr'])
			strand = pos['strand']
			tss    = pos['start'] if strand == 1 else pos['end']
			f.write ("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chr, tss, tss + 1, hit['symbol'], 0, ("+" if strand == 1 else "-"), hit['query']))
		except TypeError:
			sys.stderr.write('Encounter TypeError, hit is: %s\n' % str(hit))

