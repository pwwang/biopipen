expr    = read.table ("{{in.expr}}", sep="\t", header=TRUE, row.names = 1, check.names=F)
batch   = read.table ("{{in.batch}}", sep="\t", header=FALSE, row.names = 1, check.names=F)
batch   = factor(batch[colnames(expr),1,drop=T])
{% if args.tool | lambda x: x == 'combat' %}
library(sva)
newexpr   = ComBat(dat=expr, batch=batch, par.prior = TRUE, mod = NULL)
write.table (newexpr, "{{out.outfile}}", col.names=T, row.names=T, sep="\t", quote=F)
{% else %}
stop('Unsupported tool: {{args.tool}}.')
{% endif %}


# boxplot
{% if args.boxplot %}
{{ plotBoxplot }}
bpfile = file.path("{{out.outdir}}", "{{in.expr | fn | fn}}.boxplot.png")
plotBoxplot(newexpr, bpfile, devpars = {{args.devpars | Rlist}}, ggs = {{args.boxplotggs | Rlist}})
{% endif %}

# heatmap
{% if args.heatmap %}
{{ plotHeatmap }}
hmfile = file.path("{{out.outdir}}", "{{in.expr | fn | fn}}.heatmap.png")
plotHeatmap(newexpr, hmfile, devpars = {{args.devpars | Rlist}}, ggs = {{args.heatmapggs | Rlist}})
{% endif %}

# histgram
{% if args.histplot %}
{{ plotHist }}
histfile = file.path("{{out.outdir}}", "{{in.expr | fn | fn}}.hist.png")
plotHist(newexpr, histfile, devpars = {{args.devpars | Rlist}}, ggs = {{args.histplotggs | Rlist}})
{% endif %}
