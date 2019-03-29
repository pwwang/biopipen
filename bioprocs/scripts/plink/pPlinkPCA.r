library(methods)
{{rimport}}('__init__.r')
options(stringsAsFactors = FALSE)

indir     = {{i.indir | R}}
outfile   = {{o.outfile | R}}
outdir    = {{o.outdir | R}}
params    = {{args.params | R}}
select    = {{args.select | R}}
plots     = {{args.plots | R}}
devpars   = {{args.devpars | R}}
plink     = {{args.plink | quote}}
nthread   = {{args.nthread | R}}
samid     = {{args.samid | R}}
indep     = {{args.indep | R}}
highld    = {{args.highld | R}}
jobstdout = {{job.outfile | R}}

if (is.true(plots, 'any')) {
	{{rimport}}('plot.r')
}

bedfile = Sys.glob(file.path(indir, '*.bed'))
input   = tools::file_path_sans_ext(bedfile)
output  = file.path(outdir, basename(input))

# if (is.true(covfile) && is.true(covname)) {
# 	newcovfile = paste0(output, '.cov.txt')
# 	covdata = read.table.inopts(covfile, list(rnames = FALSE))
# 	if (covname == 'fid' || covname == 'iid') {
# 		covdata = cbind(covdata[,1,drop=FALSE], covdata)
# 	} else {
# 		ids = do.call(rbind, strsplit(as.character(covdata[,1]), '_', fixed = TRUE))
# 		nids = ncol(ids)
# 		if (nids %% 2 > 0) {
# 			stop('You have imbalanced _ in FID and IID')
# 		}
# 		midunscore = as.integer(nids/2)
# 		ids = data.frame(
# 			FID = apply(ids[, 1:midunscore], 1, function(x) paste(x, collapse = '_')),
# 			IID = apply(ids[, (midunscore + 1):nids], 1, function(x) paste(x, collapse = '_'))
# 		)
# 		covdata = cbind(ids, covdata[, -1, drop = FALSE])
# 	}
# 	colnames(covdata)[1:2] = c('FID', 'IID')
# 	write.table(covdata, newcovfile, col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
# 	covfile = newcovfile
# }

# prune high LD SNPs
pruneparams = list(
	bfile            = input,
	exclude          = highld,
	range            = T,
	`indep-pairwise` = indep,
	out              = output
)
cmd = sprintf("%s %s 1>&2", plink, cmdargs(pruneparams, equal = ' '))
runcmd(cmd)

prunein = paste0(output, '.prune.in')

default.params = list(
	pca     = TRUE,
	bfile   = input,
	out     = output,
	threads = nthread,
	extract = prunein
	#covar   = ifelse(is.true(covfile), covfile, FALSE)
)

params = update.list(default.params, params)
# stdout slows down the program
cmd = sprintf("%s %s 1>&2", plink, cmdargs(params, equal = ' '), jobstdout)
runcmd(cmd)

valfile = paste0(output, '.eigenval')
vecfile = paste0(output, '.eigenvec')

sdevdf = read.table(valfile, header = FALSE, row.names = NULL)
rownames(sdevdf) = paste0('PC', 1:nrow(sdevdf))
colnames(sdevdf) = 'Sdev'

sdevdf$Percent    = 100 * sdevdf$Sdev / sum(sdevdf$Sdev)
sdevdf$CumPercent = cumsum(sdevdf$Percent)

write.table(sdevdf, paste0(output, '.sdev.txt'), sep = "\t", quote = FALSE)

npcs = nrow(sdevdf)
if (is.null(select)) {
	select = 1:npcs
} else if (select >= 1) {
	select = 1:min(npcs, select)
} else {
	select = 1:max(2, min(npcs, sum(sdevdf$Percent >= select * 100)))
}

pcs = read.table(vecfile, header = FALSE, row.names = NULL, sep = " ", check.names = FALSE)
if (samid == 'fid') {
	rownames(pcs) = pcs[,1]
} else if (samid == 'iid') {
	rownames(pcs) = pcs[,2]
} else {
	rownames(pcs) = paste(pcs[,1], pcs[,2], sep = '_')
}
pcs = pcs[, -(1:2), drop = FALSE]
colnames(pcs) = paste0('PC', 1:ncol(pcs))
write.table(pcs, paste0(output, '.allpcs.txt'), sep = "\t", quote = FALSE)
write.table(pcs[, select, drop = FALSE], outfile, sep = "\t", quote = FALSE)

if (is.true(plots$scree)) {
	screefile = paste0(output, '.screeplot.png')
	if (!is.list(plots$scree)) {
		plots$scree = list(ncp = 20)
	}
	plot.bar(
		cbind(sdevdf[1:plots$scree$ncp,,drop = FALSE], PCS = paste0('PC', 1:plots$scree$ncp)), 
		screefile,
		x = 'PCS', y = 'Percent'
	)
}

if (is.true(plots$pairs$anno)) {
	options(stringsAsFactors = TRUE)
	annos = read.table.inopts(plots$pairs$anno, list(cnames = TRUE, rnames = TRUE))
	annos = annos[rownames(pcs),,drop = FALSE]
	pdata = cbind(pcs[, 1:plots$pairs$ncp,drop = FALSE], annos)
	pairsfile = paste0(output, '.pairs.png')
	if (is.null(plots$pairs$params))
		plots$pairs$params = list()
	if (is.null(plots$pairs$ggs))
		plots$pairs$ggs = list()
	plot.pairs(pdata, pairsfile, params = plots$pairs$params, ggs = plots$pairs$ggs)
}





