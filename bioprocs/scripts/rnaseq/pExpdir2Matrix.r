setwd("{{in.expdir}}")

isGoodRname = function(rname) {
	for (excl in {{args.exrows | Rlist}}) {
		if (grepl(excl, rname)) {
			return (FALSE)
		}
	}
	return (TRUE)
}
isGoodRname = Vectorize(isGoodRname)

exp = NULL
for (efile in Sys.glob({{args.pattern | quote}})) {
	cat("pyppl.log: Reading", efile, "...\n", file = stderr())
	sample = tools::file_path_sans_ext(basename(efile))
	if (grepl ('.gz$', efile)) efile = gzfile (efile)
	tmp    = read.table (efile, sep="\t", header={{args.header | R}}, row.names = 1, check.names=F)
	rnames = rownames(tmp)
	rnames = rnames[isGoodRname(rnames)]
	tmp    = tmp[rnames,,drop=F]
	{% if args.header | lambda x: not x %}
	colnames (tmp) = c(sample)
	{% endif %}
	exp    = if(is.null(exp)) tmp else cbind (exp, tmp)
}

write.table (exp, "{{out.outfile}}", col.names=T, row.names=T, sep="\t", quote=F)

# boxplot
{% if args.boxplot %}
{{ plotBoxplot }}
bpfile = file.path("{{out.outdir}}", "{{in.expdir | fn}}.boxplot.png")
plotBoxplot(exp, bpfile, devpars = {{args.devpars | Rlist}}, ggs = {{args.boxplotggs | Rlist}})
{% endif %}

# heatmap
{% if args.heatmap %}
{{ plotHeatmap }}
hmfile = file.path("{{out.outdir}}", "{{in.expdir | fn}}.heatmap.png")
plotHeatmap(exp, hmfile, devpars = {{args.devpars | Rlist}}, ggs = {{args.heatmapggs | Rlist}})
{% endif %}

# histgram
{% if args.histplot %}
{{ plotHist }}
histfile = file.path("{{out.outdir}}", "{{in.expdir | fn}}.hist.png")
plotHist(exp, histfile, devpars = {{args.devpars | Rlist}}, ggs = {{args.histplotggs | Rlist}})
{% endif %}

