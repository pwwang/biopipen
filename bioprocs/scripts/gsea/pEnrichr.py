from os import path
from enrichr import Enrichr
from pyppl import Box
from bioprocs.utils import alwaysList
from bioprocs.utils.gene import genenorm
from bioprocs.utils.tsvio import TsvReader

infile   = {{in.infile | quote}}
inopts   = {{args.inopts}}
genecol  = {{args.genecol | quote}}
cachedir = {{args.cachedir | quote}}
top      = {{args.top}}
title    = {{args.title | quote}}

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

en = Enrichr()
en.addList(genes, description = path.basename(infile))

dbs = alwaysList({{args.libs | quote}})
for db in dbs:
	outfile  = "{{out.outdir}}/{{in.infile | fn}}-%s.txt" % db
	en.enrich(db)
	en.export(outfile)
	{% if args.plot %}
	plotfile = "{{out.outdir}}/{{in.infile | fn}}-%s.png" % db
	en.plot(plotfile, title = title.format(library = db), top = top)
	{% endif %}
