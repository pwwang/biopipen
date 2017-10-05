
library(simpleaffy)

dir.create("{{out.outdir}}/CELs", F)
for (celfile in Sys.glob("{{in.expdir}}/{{args.pattern}}")) {
	file.link(celfile, file.path("{{out.outdir}}", "CELs", basename(celfile)))
}
setwd("{{out.outdir}}/CELs")
files    = Sys.glob("{{args.pattern}}")

gfile = "covdesc.txt"
# generate a covdesc file
if ("{{args.gfile}}" == '') {
  group  = matrix(1, ncol = 1, nrow = length(files))
  colnames(group) = c("treatment")
  rownames(group) = files
  write.table(group, gfile, quote = F, sep = "\t")
} else {
  file.link("{{args.gfile}}", gfile)
}

affydata = read.affy(covdesc = gfile)
if ("{{args.cdffile}}" != '') {
  library(makecdfenv)
  cdffile = "{{args.cdffile}}"
  cdfname = cleancdfname(whatcdf(files[1]))
  assign(cdfname, make.cdf.env(basename(cdffile), cdf.path = dirname(cdffile)))
  affydata@cdfName = cdfname
} 
exprs    = call.exprs(affydata, algorithm = "{{args.norm}}")
exprsout = exprs@assayData$exprs

# annotate
{% if args.annofile %}
annos    = read.table("{{args.annofile}}", header=F, sep="\t", row.names = 1, check.names = F )
rnames   = intersect(rownames(exprsout), rownames(annos))
exprsout = exprsout[rnames,,drop=F]
rownames(exprsout) = make.unique(as.vector(annos[rnames,,drop=T]))
{% endif %}

cnames   = colnames(exprsout)
cnames   = gsub("\\.gz$", "", cnames)
cnames   = gsub("\\.CEL$", "", cnames)
cnames   = gsub("\\.cel$", "", cnames)
colnames(exprsout) = cnames
write.table(exprsout, "{{out.outfile}}", quote=F, sep="\t")

exprsout = log2(exprsout + 1)
# boxplot
{% if args.boxplot %}
{{ plotBoxplot }}
bpfile = file.path("{{out.outdir}}", "{{in.expdir | fn}}.boxplot.png")
plotBoxplot(exprsout, bpfile, devpars = {{args.devpars | Rlist}}, ggs = {{args.boxplotggs | Rlist}})
{% endif %}

# heatmap
{% if args.heatmap %}
{{ plotHeatmap }}
hmfile = file.path("{{out.outdir}}", "{{in.expdir | fn}}.heatmap.png")
hmexp  = if (nrow(exprsout) > {{args.heatmapn}}) exprsout[sample(nrow(exprsout),size={{args.heatmapn}}),] else exprsout
plotHeatmap(hmexp, hmfile, devpars = {{args.devpars | Rlist}}, ggs = {{args.heatmapggs | Rlist}})
{% endif %}

# histgram
{% if args.histplot %}
{{ plotHist }}
histfile = file.path("{{out.outdir}}", "{{in.expdir | fn}}.hist.png")
plotHist(exprsout, histfile, devpars = {{args.devpars | Rlist}}, ggs = {{args.histplotggs | Rlist}})
{% endif %}

