{{rimport}}('__init__.r', 'plot.r')

indir     = {{i.indir | R}}
exrows    = {{args.exrows | R}}
fn2sample = Vectorize({{args.fn2sample}})
pattern   = {{args.pattern | R}}
outfile   = {{o.outfile | R}}
outdir    = {{o.outdir | R}}
plot      = {{args.plot | R}}
prefix    = {{i.indir, args.pattern | dirpat2name | R}}
devpars   = {{args.devpars | R}}
hmrows    = {{args.hmrows | R}}
ggs       = {{args.ggs | R}}

setwd(indir)

isGoodRname = function(rname) {
	for (excl in exrows) {
		if (grepl(excl, rname)) {
			return (FALSE)
		}
	}
	return (TRUE)
}
isGoodRname = Vectorize(isGoodRname)


exp = NULL
for (efile in Sys.glob(pattern)) {

	logger("Reading", efile, "...")
	
	sample = tools::file_path_sans_ext(basename(efile))
	sample = fn2sample(sample)

	if (grepl ('.gz$', efile)) efile = gzfile (efile)

	tmp = read.table.nodup (efile, sep="\t", header = F, row.names = 1, check.names = F)
	tmp = tmp[isGoodRname(rownames(tmp)), , drop=F]

	colnames(tmp) = sample

	exp    = if(is.null(exp)) tmp else cbind (exp, tmp)
}

write.table (exp, outfile, col.names = T, row.names = T, sep = "\t", quote = F)

exp = log2(exp + 1)

if (plot$boxplot) {
	bpfile = file.path(outdir, paste0(prefix, '.boxplot.png'))
	plot.boxplot(exp, bpfile, stack = T, devpars = devpars, ggs = ggs$boxplot)
}

if (plot$heatmap) {
	hmfile = file.path(outdir, paste0(prefix, ".heatmap.png"))
	hmexp  = if (nrow(exp) > hmrows) exp[sample(nrow(exp),size=hmrows),] else exp
	plot.heatmap(hmexp, hmfile, devpars = devpars, ggs = ggs$heatmap)
}

if (plot$histogram) {
	histfile = file.path(outdir, paste0(prefix, ".histo.png"))
	plot.histo(stack(as.data.frame(exp)), histfile, devpars = devpars, ggs = ggs$histogram)
}
