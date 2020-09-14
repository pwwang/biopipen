{{"__init__.R" | rimport}}

library(methods)
library(atSNP)
library(reshape2)
library(parallel)
options(stringsAsFactors = FALSE)

infile = {{i.infile | R}}
snpbed = {{i.snpbed | R}}
outfile = {{o.outfile | R}}
outdir = {{o.outdir | R}}
tflist = {{args.tflist | R}}
tfmotifs = {{args.tfmotifs | R}}
genome = {{args.genome | R}}
fdr = {{args.fdr | R}}
pval = {{args.pval | R}}
stacked = {{args.stacked | R}}
doplot = {{args.plot | R}}
nthread = {{args.nthread | R}}
devpars = {{args.devpars | R}}

genome_lib = sprintf("BSgenome.Hsapiens.UCSC.%s", genome)
library(genome_lib, character.only = TRUE)

### load desired pwms
# load the list (or motifs)
# col1: motif
# col2: tf
logger("Loading motif tf names ...")
tflist = read.table.inopts(tflist, list(cnames=FALSE, rnames=FALSE))

# load all of 'em first
logger("Loading all motifs ...")
pwms = LoadMotifLibrary(tfmotifs)
pwms = pwms[unique(tflist[, 1])]

### read snp-tf pairs
logger("Reading snp-tf pairs ...")
if (stacked) {
	snptfs = read.table.inopts(infile, list(cnames=FALSE, rnames=FALSE))
	colnames(snptfs) = c("SNP", "TF")
	snptfs = dcast(snptfs, SNP ~ TF)
	rownames(snptfs) = snptfs[, 1]
	snptfs = snptfs[, -1]
} else {
	snptfs = read.table.inopts(infile, list(cnames=TRUE, rnames=TRUE))
}
tfs = colnames(snptfs)

tfs_notinlist = setdiff(tfs, unique(unlist(tflist[, 2])))
if (length(tfs_notinlist) > 0) {
	logger("TFs not found in tflist, will be ignored:", level = "WARNING")
	logger("  ", paste(tfs_notinlist, collapse = ", "), level = "WARNING")
}
tfs = intersect(tfs, unlist(tflist[, 2]))

### load all snp coordinates
logger("Loading all snp coordinates ...")
snps = read.table.inopts(snpbed, list(cnames=FALSE, rnames=FALSE))
snps = row.apply(snps, cnames = c("snpid", "chr", "snp", "a1", "a2"), function(row) {
	ret = list()
	ret$snpid = row[, 4]
	ret$chr = ifelse(grepl("chr", row[,1]), row[, 1], paste0("chr", row[, 1]))
	ret$snp = row[, 3]
	ret$a1 = row[, 7]
	ret$a2 = unlist(strsplit(as.character(unlist(row[, 8])), ',', fixed = TRUE))[1]
	ret
})

tfdir = file.path(outdir, 'TFs')
dir.create(tfdir, showWarnings = FALSE)

do_one_tf = function(tf) {
	logger("- Handling TF: " , tf)
	mysnps = rownames(snptfs[snptfs[, tf] > 0, tf, drop = FALSE])
	mysnps_notincoord = setdiff(mysnps, snps$snpid)
	if (length(mysnps_notincoord) > 0) {
		logger("  SNPs not found in coordinate file, ignore:", level = "WARNING")
		logger("    ", paste(mysnps_notincoord, collapse = ", "), level = "WARNING")
	}
	mysnps = intersect(mysnps, snps$snpid)
	if (length(mysnps) == 0) {
		logger("  No available SNPs found, skip.", level = "WARNING")
		return (NULL)
	}
	dir.create(file.path(tfdir, tf), showWarnings = FALSE)
	motif = tflist[tflist[, 2] == tf, 1]
	snpfile = file.path(tfdir, tf, 'snpinfo.txt')
	write.table(snps[snps$snpid %in% mysnps, , drop=FALSE],
				snpfile, row.names=FALSE, quote=FALSE, sep="\t")

	snp_info = LoadSNPData(
		snpfile,
		genome.lib       = genome_lib,
		half.window.size = max(sapply(pwms, length))/4,
		default.par      = TRUE,
		mutation         = TRUE
	)

	atsnp_scores = ComputeMotifScore(pwms[motif], snp_info)
	atsnp_result = ComputePValues(motif.lib = pwms[motif],
								  snp.info = snp_info,
								  motif.scores = atsnp_scores$motif.scores)
	atsnp_result$TF = tf
	atsnp_result = atsnp_result[,
		c('snpid', 'TF', 'motif', 'pval_ref', 'pval_snp', 'pval_rank'),
		drop = FALSE]
	colnames(atsnp_result) = c('Snp', 'TF', 'Motif', 'Pval_Ref', 'Pval_Mut', 'Pval_Diff')

	if (fdr == FALSE) {
		# we can safely remove insignificant results
		atsnp_result = atsnp_result[atsnp_result$Pval_Diff < pval, , drop = FALSE]
	}

	if (doplot) {
		for (i in 1:nrow(atsnp_result)) {
			if (atsnp_result[i, 'Pval_Diff'] >= pval) next
			snp = atsnp_result[i, 'Snp']
			plotfile = file.path(tfdir, tf, paste0(snp, '-', motif, '.png'))
			match_seq <- dtMotifMatch(
				atsnp_scores$snp.tbl,
				atsnp_scores$motif.scores,
				snpids = snp,
				motifs = motif,
				motif.lib = pwms[motif])
			do.call(png, c(list(file = plotfile), devpars))
			plotMotifMatch(match_seq,  motif.lib = pwms[motif])
			dev.off()
		}
	}
	atsnp_result
}

results = do.call(rbind, mcmapply(do_one_tf, tfs, SIMPLIFY=FALSE, mc.cores=nthread))
#results = do.call(rbind, apply(tfs, do_one_tf))

if (fdr != FALSE && !is.null(results) && nrow(results) > 0) {
	results$Qval_Diff = p.adjust(results$Pval_Diff, method = ifelse(fdr == TRUE, 'BH', fdr))
	results = results[results$Pval_Diff < pval, , drop = FALSE]
}

write.table(results, outfile, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
