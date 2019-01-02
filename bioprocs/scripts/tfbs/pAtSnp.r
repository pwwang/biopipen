library(methods)
library(atSNP)
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
	# suppose no head
	# must be a bed file with 7 columns
	# bed6 + alleles (A,C,G)
	# the first is reference allele and rest are mutant alleles
	# first 2 will be used (ref and minor alleles in this case)
	snps    = read.table(snpfile, header = FALSE, row.names = NULL, sep = "\t", check.names = FALSE)
	snpfile = file.path(outdir, 'snps.txt')
	# write it to snpfile for loading
	snps2 = data.frame(snpid = snps[,4], chr = snps[,1], snp = snps[,3])
	snps2 = cbind(snps2, t(sapply(snps[,7], function(x) unlist(strsplit(as.character(x), ',', fixed = T))[1:2])))
	rm(snps)
	colnames(snps2) = c('snpid', 'chr', 'snp', 'a1', 'a2')
	write.table(snps2, file = snpfile, row.names = FALSE, quote = FALSE)
}

snp_info = LoadSNPData(
	snpfile, 
	genome.lib       = sprintf("BSgenome.Hsapiens.UCSC.%s", genome),
	half.window.size = max(sapply(pwms, length))/4,
	default.par      = TRUE,
	mutation         = TRUE
)

atsnp.scores = ComputeMotifScore(pwms, snp_info, ncores = nthread)
atsnp.result = ComputePValues(motif.lib = pwms, snp.info = snp_info, motif.scores = atsnp.scores$motif.scores, ncores = nthread)
atsnp.result = atsnp.result[order(pval_rank) & pval_rank <= pval, list(snpid, motif, pval_ref, pval_snp, pval_rank)]
colnames(atsnp.result) = c('Snp', 'Motif', 'Pval_Ref', 'Pval_Mut', 'Pval_Diff')
atsnp.result[, TF := paste(tflist[which(tflist[,1] == Motif), 2], collapse = ',')]
atsnp.result = atsnp.result[, list(Snp, TF, Motif, Pval_Ref, Pval_Mut, Pval_Diff)]
if (fdr != FALSE) {
	atsnp.result[, Qval_Diff := p.adjust(Pval_Diff, method = fdr)]
}
write.table(atsnp.result, outfile, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

if (doplot) {
	for (i in 1:nrow(atsnp.result)) {
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


