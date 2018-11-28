library(methods)
library(data.table)
{{rimport}}('__init__.r', 'plot.r')

indir  = {{i.indir | R}}
outdir = {{o.outdir | R}}

plink   = {{args.plink | R}}
indep   = {{args.indep | R}}
highld  = {{args.highld | R}}
devpars = {{args.devpars | R}}
pihat   = {{args.pihat | R}}

bedfile = Sys.glob(file.path(indir, '*.bed'))
input   = tools::file_path_sans_ext(bedfile)
output  = file.path(outdir, basename(input))

params = list(
	bfile            = input,
	exclude          = highld,
	range            = T,
	`indep-pairwise` = indep,
	out              = output
)
cmd = sprintf("%s %s 1>&2", plink, cmdargs(params, equal = ' '))
runcmd(cmd)

prunein = paste0(output, '.prune.in')
params = list(
	bfile   = input,
	extract = prunein,
	genome  = T,
	out     = output
)
cmd = sprintf("%s %s 1>&2", plink, cmdargs(params, equal = ' '))
runcmd(cmd)

genome = read.table(paste0(output, '.genome'), row.names = NULL, header = T, check.names = F)
# "unmelt" it
# FID1 IID1 FID2 IID2 RT EZ Z0     Z1     Z2     PI_HAT PHE DST      PPC    RATIO
# s1   s1   s2   s2   UN NA 1.0000 0.0000 0.0000 0.0000 -1  0.866584 0.0000 0.9194
# s1   s1   s2   s2   UN NA 0.4846 0.3724 0.1431 0.3293 -1  0.913945 0.7236 2.0375
# s1   s1   s3   s3   UN NA 1.0000 0.0000 0.0000 0.0000 -1  0.867186 0.0000 1.0791
genome$SAMPLE1 = paste(genome$FID1, genome$IID1, sep = "\t")
genome$SAMPLE2 = paste(genome$FID2, genome$IID2, sep = "\t")

# get all samples
samples    = unique(c(genome$SAMPLE1, genome$SAMPLE2))
# make paired into a distance-like matrix
similarity = dcast(genome, SAMPLE1 ~ SAMPLE2, value.var = "PI_HAT")
rm(genome)
# get the rownames back
rownames(similarity) = similarity[, 1]
similarity = similarity[, -1, drop = F]
# get samples that didn't involved
missedrow  = setdiff(samples, rownames(similarity))
missedcol  = setdiff(samples, colnames(similarity))
similarity[missedrow, ] = NA
similarity[, missedcol] = NA
# order the matrix
similarity = similarity[samples, samples, drop = F]
# transpose the matrix to get the symmetric values
sim2 = t(similarity)
isna = is.na(similarity)
# fill the na's with their symmetric values
similarity[isna] = sim2[isna]
rm(sim2)
# still missing: keep them
similarity[is.na(similarity)] = 0
# get the marks (samples that fail the pihat cutoff)
nsams = length(samples)
fails = which(similarity > pihat)
marks = data.frame(x = (fails - 1)%%nsams + 1, y = ceiling(fails/nsams))
diag(similarity) = 1

failflags = rep(F, nrow(marks))
freqs     = as.data.frame(table(factor(as.matrix(marks))))
freqs     = freqs[order(freqs$Freq, decreasing = T), 'Var1', drop = T]
ibd.fail  = c()
while (sum(failflags) < nrow(marks)) {
	samidx   = freqs[1]
	ibd.fail = c(ibd.fail, samples[samidx])
	freqs    = freqs[-1]
	sapply(1:nrow(marks), function(i) {
		if (samidx %in% marks[i,])
			failflags[i] <<- T
	})
}
writeLines(ibd.fail, con = file(paste0(output, '.ibd.fail')))

plot.heatmap(
	similarity, 
	plotfile = paste0(output, '.ibd.png'),
	devpars  = devpars,
	params   = list(dendro = 'row'),
	ggs      = list(
		theme  = list(
			axis.text.x  = element_text(angle = 90, hjust = 1),
			legend.title = element_text(size = 7),
			legend.text = element_text(size = 7),
			legend.key.size = unit(.8, 'lines'),
			legend.background = element_rect(size=0.1, linetype="solid", fill = '#EEEEFF', colour ="darkblue")
		),
		guides = list(
			fill  = guide_legend(title = "PI_HAT"),
			color = guide_legend(title = sprintf("PI_HAT >= %s", pihat), label.theme = element_blank())
		),
		geom_point = list(aes(x = x, y = y, color = "red"), data = marks, shape = 4)
	)
)