import json
import requests
import math
from hashlib import md5
genes = [line.split()[0] for line in open("{{in.infile}}") if line.strip()]
genes = list(set(genes))

if {{args.norm}}:
	genes    = sorted(genes)
	scopes   = ['symbol', 'alias']
	fields   = ['symbol']
	species  = 'human'
	uid      = md5(''.join(genes + scopes + fields + [species])).hexdigest()[:8]
	igfile   = "{{args.mgcache}}/mygeneinfo.%s" % uid
	ogfile   = "{{out.outdir}}/input.genes"

	if path.isfile(igfile):
		with open(igfile) as f, open(ogfile, 'w') as fout:
			for line in f:
				fout.write(line)
				(query, symbol)  = line.strip().split('\t')
				genes.append(symbol)
	else:
		from mygene import MyGeneInfo
		mg       = MyGeneInfo()
		mgret    = mg.getgenes (genes, scopes=scopes, fields=fields, species=species)
		with open (igfile, "w") as fout, open(ogfile, 'w') as fout2:
			for gene in mgret:
				if not 'symbol' in gene: continue
				genes.append(gene['symbol'])
				fout.write("%s\t%s\n" % (gene['query'], gene['symbol']))
				fout2.write("%s\t%s\n" % (gene['query'], gene['symbol']))
	
if {{args.plot}}:
	import matplotlib
	matplotlib.use('Agg')
	from matplotlib import pyplot as plt
	from matplotlib import gridspec
	from matplotlib import patches

## upload
ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/addList'
genes_str   = "\n".join(genes)
description = '{{in.infile | fn}}'
payload = {
    'list': (None, genes_str),
    'description': (None, description)
}

response = requests.post(ENRICHR_URL, files=payload)
if not response.ok:
    raise Exception('Error analyzing gene list')

data = json.loads(response.text)

## do enrichment
dbs = "{{args.dbs}}".split(',')
dbs = map (lambda s: s.strip(), dbs)

ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/enrich'
query_string = '?userListId=%s&backgroundType=%s'

head = ["#Rank", "Term name", "P-value", "Z-score", "Combined score", "Overlapping genes", "Adjusted p-value", "Old p-value", "Old adjusted p-value"]
topn = {{args.topn}}
for db in dbs:
	user_list_id = data['userListId']
	gene_set_library = db
	response = requests.get(
		ENRICHR_URL + query_string % (user_list_id, gene_set_library)
	)
	if not response.ok:
		raise Exception('Error fetching enrichment results against %s' % db)
	
	data    = json.loads(response.text)
	data    = data[db]
	d2plot  = []
	outfile = "{{out.outdir}}/%s.txt" % db
	fout    = open (outfile, "w")
	fout.write ("\t".join(head) + "\n")
	for i, ret in enumerate(data):
		fout.write ("\t".join(['|'.join(r) if x == 5 else str(r) for x,r in enumerate(ret)]) + "\n")
		if topn < 1 and ret[2] >= topn: continue
		if topn >= 1 and i > topn - 1: continue
		if {{args.rmtags}} and "_Homo sapiens_hsa" in ret[1]: ret[1] = ret[1][:-22]
		d2plot.append (ret)
	fout.close()
	
	if {{args.plot}}:
		#d2plot   = sorted (d2plot, cmp=lambda x,y: 0 if x[2] == y[2] else (-1 if x[2] < y[2] else 1))
		plotfile = "{{out.outdir}}/%s.png" % db
		gs = gridspec.GridSpec(1, 2, width_ratios=[3, 7]) 
		rownames = [r[1] if len(r[1])<=40 else r[1][:40] + ' ...' for r in d2plot]
		rnidx    = range (len (rownames))
		ax1 = plt.subplot(gs[0])
		plt.title ("{{args.title}}".replace("{db}", db), fontweight='bold')
		
		ax1.xaxis.grid(alpha=.6, ls = '--', zorder = -99)
		plt.subplots_adjust(wspace=.01, left=0.5)
		ax1.barh(rnidx, [len(r[5]) for r in d2plot], color='blue', alpha=.6)
		plt.yticks (rnidx, rownames)
		ax1.yaxis.set_ticks_position('none')
		ax1.tick_params(axis='x', colors='blue')
		ax1.spines['top'].set_visible(False)
		ax1.spines['left'].set_visible(False)
		ax1.spines['right'].set_visible(False)
		ax1.spines['bottom'].set_linewidth(1)
		ax1.invert_xaxis()
		ax1.invert_yaxis()
		xticks = ax1.xaxis.get_major_ticks()
		xticks[0].label1.set_visible(False)
		
		ax2 = plt.subplot(gs[1])
		ax2.xaxis.grid(alpha=.6, ls = '--', zorder = -99)
		ax2.barh(rnidx, [-math.log(r[2], 10) for r in d2plot], color='red', alpha = .6)
		for i, r in enumerate(d2plot):
			t  = str("%.2E" % r[2])
			tx = 0.1
			ty = i + 0.1
			ax2.text(tx, ty, t, fontsize=8)
		ax2.tick_params(axis='x', colors='red')
		ax2.spines['top'].set_visible(False)
		ax2.spines['left'].set_visible(False)
		ax2.spines['right'].set_visible(False)
		ax2.spines['bottom'].set_linewidth(1)
		ax2.invert_yaxis()
		plt.yticks([])
		ng_patch = patches.Patch(color='blue', alpha=.6, label='# overlapped genes')
		pv_patch = patches.Patch(color='red', alpha = .6, label='-log(p-value)')
		plt.figlegend(handles=[ng_patch, pv_patch], labels=['# overlapped genes', '-log(p-value)'], loc="lower center", ncol=2, edgecolor="none")
		plt.savefig(plotfile, dpi=300)