library(methods)
library(atSNP)
library(chunked)
options(stringsAsFactors = FALSE)

snpfile  = {{i.snpfile | R}}
tffile   = {{i.tffile | R}}
outfile  = {{o.outfile | R}}
outdir   = {{o.outdir | R}}
tfmotifs = {{args.tfmotifs | R}}
fdr      = {{args.fdr | :'BH' if a is True else a | R}}
pval     = {{args.pval | R}}
nthread  = {{args.nthread | R}}
genome   = {{args.genome | R}}
doplot   = {{args.plot | R}}
devpars  = {{args.devpars | R}}
prefix   = file.path(outdir, {{o.outfile | fn2 | R}})

logger = function(..., level = 'LOG') {
	msg = paste(...)
	if (!endsWith(msg, "\n")) msg = paste0(msg, "\n")
	cat(paste0('pyppl.log.', level, ': ', msg), file = stderr())
}

### load desired pwms
# load the list (or motifs)
# col1: motif
# col2: tf
tflist = read.table(tffile, header = FALSE, row.names = NULL, sep = "\t", check.names = FALSE)
# load all of 'em first
pwms = LoadMotifLibrary(tfmotifs)
pwms = pwms[unique(tflist[,1])]

### load snp information
if (endsWith(snpfile, '.bed')) {
	logger('Reformatting bed file to atSNP format ...')
	# suppose no head
	# must be a bed file with 8 columns
	# bed6 + ref + alt alleles (C,G)
	# first 2 will be used (ref and minor alleles in this case)
	snps    = read.table(snpfile, header = FALSE, row.names = NULL, sep = "\t", check.names = FALSE)
	snpfile = file.path(outdir, 'snps.txt')
	# write it to snpfile for loading
	snps2 = data.frame(snpid = snps[,4], chr = snps[,1], snp = snps[,3], a1 = snps[,7])
	snps2 = cbind(snps2, a2 = sapply(snps[,8], function(x) unlist(strsplit(as.character(x), ',', fixed = TRUE))[1]))
	rm(snps)
	write.table(snps2, file = snpfile, row.names = FALSE, quote = FALSE, sep = "\t")
}

chunksize = 50
logger('Handling the snp file in chunks (chunksize:', chunksize, ')')
chunks  = read_chunkwise(snpfile, chunk_size = chunksize, format = "table", header = TRUE)
i       = 0
results = NULL
while (!chunks$is_complete()) {
	i = i + 1
	logger('- Handling chunk #', i, '...')
	sfile = file.path(outdir, paste0('snp-chunk', i, '.txt'))
	nextChunk = chunks$next_chunk()
	if (is.null(nextChunk)) {
		break
	}
	write.table(nextChunk, sfile, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
	snp_info = LoadSNPData(
		sfile, 
		genome.lib       = sprintf("BSgenome.Hsapiens.UCSC.%s", genome),
		half.window.size = max(sapply(pwms, length))/4,
		default.par      = TRUE,
		mutation         = TRUE
	)
	file.remove(sfile)

	atsnp.scores = ComputeMotifScore(pwms, snp_info, ncores = nthread)
	atsnp.result = ComputePValues(motif.lib = pwms, snp.info = snp_info, motif.scores = atsnp.scores$motif.scores, ncores = nthread)
	if (fdr == FALSE) {
		atsnp.result = atsnp.result[order(pval_rank) & pval_rank < pval, list(snpid, motif, pval_ref, pval_snp, pval_rank)]
	}
	colnames(atsnp.result) = c('Snp', 'Motif', 'Pval_Ref', 'Pval_Mut', 'Pval_Diff')
	atsnp.result[, TF := paste(tflist[which(tflist[,1] == Motif), 2], collapse = ','), by = Motif]
	atsnp.result = atsnp.result[, list(Snp, TF, Motif, Pval_Ref, Pval_Mut, Pval_Diff)]
	if (is.null(results)) {
		results = atsnp.result
	} else {
		results = rbind(results, atsnp.result)
	}

	if (doplot) {
		for (i in 1:nrow(atsnp.result)) {
			if (atsnp.result[i, Pval_Diff] >= pval) next
			plotfile = paste(prefix, atsnp.result[i, Snp], atsnp.result[i, TF], atsnp.result[i, Motif], 'png', sep = '.')
			devpars$file = plotfile
			do.call(png, devpars)
			plotMotifMatch(
				snp.tbl      = atsnp.scores$snp.tbl,
				motif.scores = atsnp.scores$motif.scores,
				snpid        = atsnp.result[i, Snp],
				motif.lib    = pwms,
				motif        = atsnp.result[i, Motif])
			dev.off()
		}
	}
}
if (fdr != FALSE) {
	results[, Qval_Diff := p.adjust(Pval_Diff, method = ifelse(fdr == TRUE, 'BH', fdr))]
	results = results[order(Pval_Diff) & Pval_Diff < pval]
}

write.table(results, outfile, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)



