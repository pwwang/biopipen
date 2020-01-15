library(methods)
library(data.table)
{{rimport}}('__init__.r', 'plot.r')

indir  = {{i.indir | R}}
outdir = {{o.outdir | R}}

plink    = {{args.plink | R}}
indep    = {{args.indep | R}}
highld   = {{args.highld | R}}
devpars  = {{args.devpars | R}}
pihat    = {{args.pihat | R}}
samid    = {{R(args.samid) if not args.samid.startswith("function") else args.samid}}
annofile = {{args.anno | R}}
doplot   = {{args.plot | R}}
seed     = {{args.seed | R}}

shell$load_config(plink = plink)
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
shell$plink(params, .raise = TRUE, .report = TRUE, .fg = TRUE)$reset()

prunein = paste0(output, '.prune.in')
params = list(
	bfile   = input,
	extract = prunein,
	genome  = T,
	out     = output
)
shell$plink(params, .raise = TRUE, .report = TRUE, .fg = TRUE)$reset()

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
samids = unlist(similarity[, 1, drop = TRUE])
rownames(similarity) = samids
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
if (!is.null(ibd.fail))
	writeLines(ibd.fail, con = file(paste0(output, '.ibd.fail')))

if (doplot) {
	set.seed(seed)
	library(ComplexHeatmap)
	fontsize8 = gpar(fontsize = 8)
	fontsize9 = gpar(fontsize = 9)
	ht_opt$heatmap_row_names_gp    = fontsize8
	ht_opt$heatmap_column_names_gp = fontsize8
	ht_opt$legend_title_gp         = fontsize9
	ht_opt$legend_labels_gp        = fontsize8
	ht_opt$simple_anno_size        = unit(3, "mm")

	samids = sapply(samples, function(sid) {
		fidiid = unlist(strsplit(sid, "\t", fixed = TRUE))
		if (is.function(samid)) {
			do.call(samid, as.list(fidiid))
		} else {
			gsub("fid", fidiid[1], gsub("iid", fidiid[2], samid))
		}
	})
	rownames(similarity) = samids
	colnames(similarity) = samids

	annos = list()
	if (is.true(annofile)) {
		options(stringsAsFactors = TRUE)
		andata = read.table.inopts(annofile, list(cnames = TRUE, rnames = TRUE))
		andata = andata[samids,,drop = FALSE]
		for (anname in colnames(andata)) {
			annos[[anname]] = as.matrix(andata[, anname])
		}
		annos$annotation_name_gp = fontsize8
		annos = do.call(HeatmapAnnotation, annos)
	}

	params = list(
		name = "PI_HAT",
		cell_fun = function(j, i, x, y, width, height, fill) {
			if (similarity[i, j] > pihat && i != j)
				grid.points(x, y, pch = 4, size = unit(.5, "char"))
		},
		#heatmap_legend_param = list(
		#	title_gp  = fontsize9,
		#	labels_gp = fontsize8
		#),
		clustering_distance_rows = function(m) as.dist(1-m),
		clustering_distance_columns = function(m) as.dist(1-m),
		top_annotation = if (length(annos) == 0) NULL else annos
	)

	plot.heatmap2(
		similarity,
		paste0(output, '.ibd.png'),
		params = params,
		draw = list(
			annotation_legend_list = list(
				Legend(
					labels = paste(">", pihat),
					title = "",
					type = "points",
					pch = 4,
					title_gp = fontsize9,
					labels_gp = fontsize8)),
			merge_legend = TRUE),
		devpars = devpars)
}
