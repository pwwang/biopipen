{{rimport}}('plot.r')

library(methods)
library(simpleaffy)

indir    = {{in.indir | R}}
pattern  = {{args.pattern | R}}
norm     = {{args.norm | R}}
gfile    = {{args.gfile | R}}
cdffile  = {{args.cdffile | R}}
annofile = {{args.annofile | R}}
hmrows    = {{args.hmrows | R}}
plot     = {{args.plot | R}}
ggs      = {{args.ggs | R}}
devpars  = {{args.devpars | R}}
outdir   = {{out.outdir | R}}
outfile  = {{out.outfile | R}}
celdir   = file.path(outdir, 'CELs')
prefix   = {{in.indir, args.pattern | dirpat2name | R}}

dir.create(celdir, F)
for (celfile in Sys.glob(file.path(indir, pattern))) {
	file.link(celfile, file.path(celdir, basename(celfile)))
}
setwd(celdir)
files    = Sys.glob(pattern)

if (gfile == '') {
	# generate a covdesc file
	gfile = "covdesc.txt"
	group  = matrix(1, ncol = 1, nrow = length(files))
	colnames(group) = c("treatment")
	rownames(group) = files
	write.table(group, gfile, quote = F, sep = "\t")
} else {
	file.link(gfile, 'covdesc.txt')
	gfile = 'covdesc.txt'
}

affydata = read.affy(covdesc = gfile)
if (cdffile != '') {
  library(makecdfenv)
  cdfname = cleancdfname(whatcdf(files[1]))
  assign(cdfname, make.cdf.env(basename(cdffile), cdf.path = dirname(cdffile)))
  affydata@cdfName = cdfname
}
exprs    = call.exprs(affydata, algorithm = norm)
exprsout = exprs@assayData$exprs

# annotate
if (annofile != '') {
	annos    = read.table(annofile, header=F, sep="\t", row.names = 1, check.names = F )
	rnames   = intersect(rownames(exprsout), rownames(annos))
	exprsout = exprsout[rnames,,drop=F]
	rownames(exprsout) = make.unique(as.vector(annos[rnames,,drop=T]))
}

cnames   = colnames(exprsout)
cnames   = gsub("\\.gz$", "", cnames)
cnames   = gsub("\\.CEL$", "", cnames, ignore.case = T)
colnames(exprsout) = cnames
write.table(round(exprsout, 3), outfile, quote=F, sep="\t")

exprsout = log2(exprsout + 1)
# boxplot
if (plot$boxplot) {
	bpfile = file.path(outdir, paste0(prefix, ".boxplot.png"))
	plot.boxplot(exprsout, bpfile, stack = T, devpars = devpars, ggs = ggs$boxplot)
}

# heatmap
if (plot$heatmap) {
	hmfile = file.path(outdir, paste0(prefix, ".heatmap.png"))
	hmexp  = if (nrow(exprsout) > hmrows) exprsout[sample(nrow(exprsout),size=hmrows),] else exprsout
	plot.heatmap(hmexp, hmfile, devpars = devpars, ggs = ggs$heatmap)
}

# histgram
if (plot$histogram) {
	histfile = file.path(outdir, paste0(prefix, ".histo.png"))
	plot.histo(stack(as.data.frame(exprsout)), histfile, devpars = devpars, ggs = ggs$histogram)
}
