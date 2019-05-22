library(methods)
{{rimport}}('__init__.r', 'plot.r')

indir   = {{i.indir | R}}
outdir  = {{o.outdir | R}}
params  = {{args.params | R}}
cutoff  = {{args.cutoff | R}}
plots   = {{args.plot | R}}
devpars = {{args.devpars | R}}
plink   = {{args.plink | quote}}
nthread = {{args.nthread | R}}
bedfile = Sys.glob(file.path(indir, '*.bed'))
input   = tools::file_path_sans_ext(bedfile)
output  = file.path(outdir, basename(input))

shell$load_config(plink = plink)

params$bfile   = input
params$out     = output
params$threads = nthread

shell$plink(params, .raise = TRUE, .report = TRUE, .fg = TRUE)$reset()

sexcheck_result = paste0(output, '.sexcheck')
if (file.exists(sexcheck_result)) {
	sexcheck = read.table(sexcheck_result, header = T, row.names = NULL, check.names = F)
	sex.sample.fail = sexcheck[which(sexcheck$STATUS == 'PROBLEM'), c('FID', 'IID'), drop=F]
	write.table(sex.sample.fail, paste0(output, '.sex.fail'), col.names = F, row.names = F, sep = "\t", quote = F)
}

hwe_result = paste0(output, '.hwe')
if (file.exists(hwe_result)) {
	hardy = read.table(paste0(output, '.hwe'), header = T, row.names = NULL, check.names = F)
	if (!is.null(cutoff$hardy.hwe)) {
		hardy.fail = hardy[which(hardy$P < cutoff$hardy.hwe), 'SNP', drop = F]
		write.table(hardy.fail, paste0(output, '.hardy.fail'), col.names = F, row.names = F, sep = "\t", quote = F)
	}
	if (is.true(plots$hardy.hwe)) {
		hardy$Pval   = -log10(hardy$P)
		hardy$Status = "Pass"
		ggs          = list()
		if (!is.null(cutoff$hardy.hwe)) {
			hardy[which(hardy$SNP %in% hardy.fail$SNP), "Status"] = "Fail"
			ggs$geom_vline = list(xintercept = -log10(cutoff$hardy.hwe), color = "red", linetype="dashed")
			ggs$geom_text  = list(
				aes(x = -log10(cutoff$hardy.hwe), y = Inf, label = cutoff$hardy.hwe),
				colour="red", angle=90, vjust = 1.2, hjust = 1.2
			)
		}
		ggs$xlab  = list("-log10(HWE p-value)")
		ggs$ylab  = list("Count")
		ggs$theme = list(legend.position = "none")
		plot.histo(
			data     = hardy,
			x        = 'Pval',
			plotfile = paste0(output, '.hardy.png'),
			params   = list(aes(fill=Status), bins = 50),
			ggs      = ggs,
			devpars  = devpars
		)
	}
	if (!is.null(cutoff$hardy.mingt) || is.true(plots$hardy.mingt)) {
		mingt = data.frame(SNP = hardy$SNP)
		gts   = t(as.data.frame(lapply(hardy$GENO, function(x) as.numeric(unlist(strsplit(as.character(x), "/", fixed = TRUE))))))
		mingt$MINGT_NO = apply(gts, 1, function(x) min(x))
		mingt$MINGT    = apply(gts, 1, function(x) min(x)/sum(x))

		if (!is.null(cutoff$hardy.mingt)) {
			if (cutoff$hardy.mingt < 1) {
				mingt.fail = mingt[which(mingt$MINGT < cutoff$hardy.mingt), 'SNP', drop = F]
				write.table(mingt.fail, paste0(output, '.mingt.fail'), col.names = F, row.names = F, sep = "\t", quote = F)
			} else {
				mingt.fail = mingt[which(mingt$MINGT_NO < cutoff$hardy.mingt), 'SNP', drop = F]
				write.table(mingt.fail, paste0(output, '.mingt.fail'), col.names = F, row.names = F, sep = "\t", quote = F)
			}
		}

		if (is.true(plots$hardy.mingt)) {
			mingt$Status = "Pass"
			ggs = list()
			if (!is.null(cutoff$hardy.mingt)) {
				mingt[which(mingt$SNP %in% mingt.fail$SNP), "Status"] = "Fail"
				ggs$geom_vline = list(xintercept = cutoff$hardy.mingt, color = "red", linetype="dashed")
				ggs$geom_text  = list(
					aes(x = cutoff$hardy.mingt, y = Inf, label = cutoff$hardy.mingt),
					colour="red", angle=90, vjust = 1.2, hjust = 1.2
				)
			}
			ggs$xlab  = list("Min_GT")
			ggs$ylab  = list("Count")
			ggs$theme = list(legend.position = "none")
			plot.histo(
				data     = mingt,
				x        = 'MINGT',
				plotfile = paste0(output, '.mingt_rate.png'),
				params   = list(aes(fill=Status), bins = 50),
				ggs      = ggs,
				devpars  = devpars
			)
			plot.histo(
				data     = mingt,
				x        = 'MINGT_NO',
				plotfile = paste0(output, '.mingt_no.png'),
				params   = list(aes(fill=Status), bins = 50),
				ggs      = ggs,
				devpars  = devpars
			)
		}
	}
}

het_result = paste0(output, '.het')
if (file.exists(het_result)) {
	phet = read.table(het_result, header = T, row.names = NULL, check.names = F)
	het = data.frame(Het = 1 - phet[, "O(HOM)"]/phet[, "N(NM)"])
	rownames(het) = paste(phet$FID, phet$IID, sep = "\t")
	if (!is.null(cutoff$het)) {
		het.mean = mean(het$Het, na.rm = T)
		het.sd   = sd(het$Het, na.rm = T)
		het.fail = rownames(het[!is.na(het$Het) & (het$Het < het.mean-cutoff$het*het.sd | het$Het > het.mean+cutoff$het*het.sd),, drop = F])
		writeLines(het.fail, con = file(paste0(output, '.het.fail')))
	}
	if (is.true(plots$het)) {
		het$Status = "Pass"
		ggs        = list()
		if (!is.null(cutoff$het)) {
			het[het.fail, "Status"] = "Fail"
			ggs$geom_vline = list(xintercept = c(het.mean-cutoff$het*het.sd, het.mean+cutoff$het*het.sd), color = "red", linetype="dashed")
			ggs$geom_text  = list(
				aes(x = het.mean-cutoff$het*het.sd, y = Inf, label = sprintf('mean - %ssd (%.3f)', cutoff$het, het.mean - cutoff$het*het.sd)),
				colour="red", angle=90, vjust = 1.2, hjust = 1.2
			)
			ggs$geom_text  = list(
				aes(x = het.mean+cutoff$het*het.sd, y = Inf, label = sprintf('mean + %ssd (%.3f)', cutoff$het, het.mean + cutoff$het*het.sd)),
				colour="red", angle=90, vjust = 1.2, hjust = 1.2
			)
		}
		ggs$xlab       = list("Sample Heterozygosity")
		ggs$ylab       = list("Count")
		ggs$geom_vline = list(xintercept = het.mean, color = "blue", linetype="dashed")
		ggs$theme      = list(legend.position = "none")
		ggs$geom_text  = list(
			aes(x = het.mean, y = Inf, label = sprintf('mean (%.3f)', het.mean)),
			colour="blue", vjust = 1.5, hjust = -.1
		)
		plot.histo(
			data     = het,
			plotfile = paste0(output, '.het.png'),
			params   = list(aes(fill=Status), bins = 50),
			ggs      = ggs,
			devpars  = devpars
		)
	}
}

freq_result = paste0(output, '.frq')
if (file.exists(freq_result)) {
	freq = read.table(freq_result, header = T, row.names = NULL, check.names = F)
	if (!is.null(cutoff$freq)) {
		freq.fail = freq[which(freq$MAF < cutoff$freq), 'SNP', drop = F]
		write.table(freq.fail, paste0(output, '.freq.fail'), col.names = F, row.names = F, sep = "\t", quote = F)
	}
	if (is.true(plots$freq)) {
		freq$Status = "Pass"
		ggs         = list()
		if (!is.null(cutoff$freq)) {
			freq[which(freq$SNP %in% freq.fail$SNP), "Status"] = "Fail"
			ggs$geom_vline = list(xintercept = cutoff$freq, color = "red", linetype="dashed")
			ggs$geom_text  = list(
				aes(x = cutoff$freq, y = Inf, label = cutoff$freq),
				colour="red", angle=90, vjust = 1.2, hjust = 1.2
			)
		}
		ggs$xlab       = list("MAF")
		ggs$ylab       = list("Count")
		ggs$theme      = list(legend.position = "none")
		plot.histo(
			data     = freq,
			x        = 'MAF',
			plotfile = paste0(output, '.freq.png'),
			params   = list(aes(fill=Status), bins = 50),
			ggs      = ggs,
			devpars  = devpars
		)
	}
}

imiss_result = paste0(output, '.imiss')
if (file.exists(imiss_result)) {
	imiss = read.table(imiss_result, header = T, row.names = NULL, check.names = F)
	callrate.sample = data.frame(Callrate = 1-imiss$F_MISS)
	rownames(callrate.sample) = paste(imiss$FID, imiss$IID, sep = "\t")
	if (!is.null(cutoff$missing.sample)) {
		callrate.sample.fail = rownames(callrate.sample[callrate.sample$Callrate < cutoff$missing.sample, , drop = F])
		writeLines(callrate.sample.fail, con = file(paste0(output, '.samplecr.fail')))
	}
	if (is.true(plots$missing.sample)) {
		callrate.sample$Status = "Pass"
		ggs = list()
		if (!is.null(cutoff$missing.sample)) {
			callrate.sample[callrate.sample.fail, "Status"] = "Fail"
			ggs$geom_vline = list(xintercept = cutoff$missing.sample, color = "red", linetype="dashed")
			ggs$geom_text  = list(
				aes(x = cutoff$missing.sample, y = Inf, label = cutoff$missing.sample),
				colour="red", angle=90, vjust = 1.2, hjust = 1.2
			)
		}
		ggs$xlab  = list("Sample Call Rate")
		ggs$ylab  = list("Count")
		ggs$theme = list(legend.position = "none")
		plot.histo(
			data     = callrate.sample,
			plotfile = paste0(output, '.samplecr.png'),
			params   = list(aes(fill=Status), bins = 50),
			ggs      = ggs
		)
	}
}

lmiss_result = paste0(output, '.lmiss')
if (file.exists(lmiss_result)) {
	lmiss = read.table(lmiss_result, header = T, row.names = NULL, check.names = F)
	lmiss$Callrate = 1-lmiss$F_MISS
	if (!is.null(cutoff$missing.snp)) {
		callrate.snp.fail = lmiss[which(lmiss$Callrate < cutoff$missing.snp), 'SNP', drop = F]
		write.table(callrate.snp.fail, paste0(output, '.snpcr.fail'), row.names = F, col.names = F, sep = "\t", quote = F)
	}
	if (is.true(plots$missing.snp)) {
		lmiss$Status = "Pass"
		ggs = list()
		if (!is.null(cutoff$missing.snp)) {
			lmiss[which(lmiss$Callrate < cutoff$missing.snp), "Status"] = "Fail"
			ggs$geom_vline = list(xintercept = cutoff$missing.snp, color = "red", linetype="dashed")
			ggs$geom_text  = list(
				aes(x = cutoff$missing.snp, y = Inf, label = cutoff$missing.snp),
				colour="red", angle=90, vjust = 1.2, hjust = 1.2
			)
		}
		ggs$xlab  = list("SNP Call Rate")
		ggs$ylab  = list("Count")
		ggs$theme = list(legend.position = "none")
		plot.histo(
			data     = lmiss,
			plotfile = paste0(output, '.snpcr.png'),
			x        = 'Callrate',
			params   = list(aes(fill=Status), bins = 50),
			ggs      = ggs,
			devpars  = devpars
		)
	}
}

