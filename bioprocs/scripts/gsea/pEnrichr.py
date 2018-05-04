import json
import requests
import math
from hashlib import md5
from pyppl import Box
from bioprocs.utils import alwaysList
from bioprocs.utils.gene import genenorm
from bioprocs.utils.tsvio import TsvReader, TsvWriter, TsvRecord
{% if args.plot %}
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt
from matplotlib import gridspec
from matplotlib import patches
{% endif %}

infile   = {{in.infile | quote}}
inopts   = {{args.inopts}}
genecol  = {{args.genecol | quote}}
cachedir = {{args.cachedir | quote}}
topn     = {{args.topn}}

{% if args.norm %}

gmapfile = "{{out.outdir}}/{{in.infile | bn}}.gnorm"
gmap = genenorm(
	infile   = infile,
	outfile  = gmapfile,
	inopts   = inopts,
	outopts  = {'head': False},
	genecol  = genecol,
	cachedir = cachedir
)
genes = [g.symbol for g in gmap.values()]

{% else %}

reader = TsvReader(infile, **inopts)
genes  = [r[genecol] for r in reader]

{% endif %}

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
	import re
	msg = re.sub('<[^<]+?>', '', response.text).splitlines()
	msg = [line for line in msg if line.strip()]
	msg = '\n'.join(msg)
	raise Exception('Error analyzing gene list: %s' % msg)

data = json.loads(response.text)

## do enrichment
dbs = alwaysList({{args.dbs | quote}})

ENRICHR_URL = 'http://amp.pharm.mssm.edu/Enrichr/enrich'
query_string = '?userListId=%s&backgroundType=%s'

head = ["Rank", "Term", "Pval", "Zscore", "CombinedScore", "OverlappingGenes", "AdjustedPval", "OldPval", "OldAdjustedPval"]

for db in dbs:
	user_list_id = data['userListId']
	gene_set_library = db
	response = requests.get(
		ENRICHR_URL + query_string % (user_list_id, gene_set_library)
	)
	if not response.ok:
		import re
		msg = re.sub('<[^<]+?>', '', response.text).splitlines()
		msg = [line for line in msg if line.strip()]
		msg = '\n'.join(msg)
		raise Exception('Error fetching enrichment results against %s: %s' % (db, msg))

	data    = json.loads(response.text)
	data    = data[db]
	d2plot  = []
	outfile = "{{out.outdir}}/{{in.infile | fn}}-%s.txt" % db
	writer  = TsvWriter(outfile)
	writer.meta.add(*head)
	writer.writeHead()
	for i, ret in enumerate(data):
		r = TsvRecord()
		r.Rank             = ret[0]
		r.Term             = ret[1]
		r.Pval             = ret[2]
		r.Zscore           = ret[3]
		r.CombinedScore    = ret[4]
		r.OverlappingGenes = '|'.join(ret[5])
		r.AdjustedPval     = ret[6]
		r.OldPval          = ret[7]
		r.OldAdjustedPval  = ret[8]
		writer.write(r)
		if topn < 1 and ret[2] >= topn: continue
		if topn >= 1 and i > topn - 1: continue
		if {{args.rmtags}} and "_" in ret[1]: ret[1] = ret[1].split('_')[0]
		d2plot.append (ret)
	writer.close()

	if {{args.plot}}:
		#d2plot   = sorted (d2plot, cmp=lambda x,y: 0 if x[2] == y[2] else (-1 if x[2] < y[2] else 1))
		plotfile = "{{out.outdir}}/{{in.infile | fn}}-%s.png" % db
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
