{{rimport}}('plot.r')

infile  = {{in.infile | R}}
outfile = {{out.outfile | R}}
inunit  = {{args.inunit | R}}
outunit = {{args.outunit | R}}
meanfl  = {{args.meanfl | R}}
nreads  = {{args.nreads | R}}
refexon = {{args.refexon | R}}
inform  = {{args.inform | lambda x: 'NULL' if not x else x}}
outform = {{args.outform | lambda x: 'NULL' if not x else x}}

data    = read.table.nodup (infile, sep="\t", header=T, row.names=1, check.names=F)
if (!is.null(inform)) {
	data = inform(data)
}

glenFromExon = function(exonfile) {
	gff  = read.table(exonfile, header = F, row.names = NULL)
	# V4: start, V5: end, V10: gene name
	glen = aggregate(V5-V4+1 ~ V10, gff, sum)
	genes = glen[,1]
	glen = glen[,-1,drop=F]
	rownames(glen) = genes
	glen
}

meanflFromFile = function(samples, mflfile) {
	if (is.numeric(mflfile)) {
		ret = matrix(mflfile, nrow = length(samples), ncol = 1)
		rownames(ret) = samples
	} else {
		ret = read.table(mflfile, header = F, row.names = 1, check.names = F, sep = "\t")
		ret = ret[samples,,drop = F]
	}
	ret
}

nreadsFromFile = function(samples, nreads) {
	if (is.numeric(nreads)) {
		ret = matrix(nreads, nrow = length(samples), ncol = 1)
		rownames(ret) = samples
	} else {
		ret = read.table(nreads, header = F, row.names = 1, check.names = F, sep = "\t")
		ret = ret[samples,,drop = F]
	}
	ret
}

# raw counts to something else
if (inunit == 'count' || inunit == 'rawcount' || inunit == 'rawcounts' || inunit == 'counts') {

	if (outunit == 'cpm') {
		samples = colnames(data)
		exp     = sapply(samples, function(s){
			# log(cpm) = log(counts) - log(sum(counts)) + log(1e6)
			# cpm = exp( log(counts) - log(sum(counts)) + log(1e6) )
			N = sum(data[, s])
			exp( log(data[, s]) - log(sum(data[, s])) + log(1e6) )
		})
		rownames(exp) = rownames(data)

	} else if (outunit == 'fpkm' || outunit == 'rpkm') {
		
		glen     = glenFromExon(refexon)
		refgenes = rownames(glen)
		mygenes  = rownames(data)
		outgenes = intersect(refgenes, mygenes)
		if (length(outgenes) < length(mygenes))
			logger('Genes not found in refexon: ', paste(setdiff(mygenes, outgenes), collapse = '\n'), level = 'WARNING')
		data    = data[outgenes, , drop = F]
		glen    = glen[outgenes, , drop = T]

		samples = colnames(data)
		fld     = meanflFromFile(samples, meanfl)

		exp     = sapply(samples, function(s){
			# N <- sum(counts)
			# exp( log(counts) + log(1e9) - log(effLen) - log(N) )
			N = sum(data[, s])
			exp( log(data[, s]) + log(1e9) - log(glen - fld[s, ] + 1) - log(N) )
		})
		rownames(exp) = outgenes

	} else if (outunit == 'fpkm-uq' || outunit == 'rpkm-uq') {
		glen     = glenFromExon(refexon)
		refgenes = rownames(glen)
		mygenes  = rownames(data)
		outgenes = intersect(refgenes, mygenes)
		if (length(outgenes) < length(mygenes))
			logger('Genes not found in refexon: ', paste(setdiff(mygenes, outgenes), collapse = '\n'), level = 'WARNING')
		data    = data[outgenes, , drop = F]
		glen    = glen[outgenes, , drop = T]

		samples = colnames(data)
		fld     = meanflFromFile(samples, meanfl)

		exp     = sapply(samples, function(s){
			# N <- sum(counts)
			# exp( log(counts) + log(1e9) - log(effLen) - log(N) )
			# N = sum(data[, s])
			RC75 = quantile(data[, s], .75)
			exp( log(data[, s]) + log(1e9) - log(glen - fld[s, ] + 1) - log(RC75) )
		})
		rownames(exp) = outgenes

	} else if (outunit == 'tpm') {
		glen     = glenFromExon(refexon)
		refgenes = rownames(glen)
		mygenes  = rownames(data)
		outgenes = intersect(refgenes, mygenes)
		if (length(outgenes) < length(mygenes))
			logger('Genes not found in refexon: \n', paste(setdiff(mygenes, outgenes), collapse = '\n'), level = 'WARNING')
		data    = data[outgenes, , drop = F]
		glen    = glen[outgenes, , drop = T]

		samples = colnames(data)
		fld     = meanflFromFile(samples, meanfl)

		exp     = sapply(samples, function(s){
			# rate <- log(counts) - log(effLen)
			# denom <- log(sum(exp(rate)))
			# exp(rate - denom + log(1e6))
			rate  = log(data[, s]) - log(glen - fld[s, ] + 1)
			denom = log(sum(exp(rate)))
			exp(rate - denom + log(1e6))
		})
		rownames(exp) = outgenes

	} else if (outunit == 'tmm') {
		library('coseq')
		exp = transform_RNAseq(data, norm="TMM")
		exp = exp$normCounts
	} else {
		stop(paste("Don't know yet how to convert raw counts to", outunit))
	}

} else if (inunit == 'fpkm' || inunit == 'rpkm') {

	if (outunit == 'count' || outunit == 'rawcount' || outunit == 'rawcounts' || outunit == 'counts') {
		glen     = glenFromExon(refexon)
		refgenes = rownames(glen)
		mygenes  = rownames(data)
		outgenes = intersect(refgenes, mygenes)
		if (length(outgenes) < length(mygenes))
			logger('Genes not found in refexon: ', paste(setdiff(mygenes, outgenes), collapse = '\n'), level = 'WARNING')
		data    = data[outgenes, , drop = F]
		glen    = glen[outgenes, , drop = T]

		samples = colnames(data)
		fld     = meanflFromFile(samples, meanfl)
		totalnr = nreadsFromFile(samples, nreads)

		exp     = sapply(samples, function(s){
			# raw count to fpkm:
			# N <- nreads
			# fpkm = exp( log(counts) + log(1e9) - log(effLen) - log(N) )
			# log(counts) = log(fpkm) + log(N) + log(effLen) - log(1e9)
			# counts = exp( log(fpkm) + log(N) + log(effLen) - log(1e9) )
			N = totalnr[s, ]
			exp( log(data[, s]) + log(N) + log(glen - fld[s, ] + 1) - log(1e9) )
		})
		rownames(exp) = outgenes
	} else if (outunit == 'tpm') {
		# tpm = exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
		exp = sapply(samples, function(s) {
			exp( log(data[, s]) - log(sum(data[, s])) + log(1e6) )
		})
		rownames(exp) = rownames(data)
	} else if (outunit == 'cpm') {
		# log(cpm) = log(counts) - log(N) + log(1e6)
		# log(fpkm) = log(counts) + log(1e9) - log(effLen) - log(N) 
		# log(cpm) - log(fpkm) = log(1e6) - log(1e9) - log(effLen)
		# cpm = exp( log(fpkm) - log(1e3) - log(effLen) )

		glen     = glenFromExon(refexon)
		refgenes = rownames(glen)
		mygenes  = rownames(data)
		outgenes = intersect(refgenes, mygenes)
		if (length(outgenes) < length(mygenes))
			logger('Genes not found in refexon: ', paste(setdiff(mygenes, outgenes), collapse = '\n'), level = 'WARNING')
		data    = data[outgenes, , drop = F]
		glen    = glen[outgenes, , drop = T]

		samples = colnames(data)

		exp = sapply(samples, function(s) {
			exp( log(data[, s]) - log(1e3) - log(glen - fld[s, ] + 1) )
		})
		rownames(exp) = outgenes

	} else {
		stop(paste("Don't know yet how to convert fpkm/rpkm to", outunit))
	}

} else if (inunit == 'tpm') {

	if (outunit == 'count' || outunit == 'rawcount' || outunit == 'rawcounts' || outunit == 'counts') {

		samples = colnames(data)
		totalnr = nreadsFromFile(samples, nreads)
		ngenes  = nrow(data)

		exp     = sapply(samples, function(s){
			# counts to tpm:
			# rate <- log(counts) - log(effLen)
			# denom <- log(sum(exp(rate)))
			# tpm = exp(rate - denom + log(1e6))
			# so:
			# log(tpm) = rate - denom + log(1e6)
			# rate = log(tpm) + denom - log(1e6)
			# log(counts) - log(effLen) = log(tpm) + log(sum(exp(rate))) - log(1e6)
			# log(counts) - log(effLen) = log(tpm) + log(sum(exp(log(counts) - log(effLen)))) - log(1e6)
			# log(counts) - log(effLen) = log(tpm) + log(sum(exp(log(counts))/exp(log(effLen)))) - log(1e6)
			# log(counts) - log(effLen) = log(tpm) + log(sum(counts/effLen)) - log(1e6)
			#                                                ?????????????
			# ??? estimated by sum(counts)/sum(effLen) * length(effLen)
			# log(counts) = log(effLen) + log(tpm) + log(sum(counts)) - log(effLen) + log(length(effLen))) - log(1e6)
			# counts = expr( log(tpm) + log(nreads) + log(length(effLen)) - log(1e6) )
			exp( log(data[, s]) + log(totalnr[s, ]) + log(ngenes) - log(1e6) )
		})
		rownames(exp) = rownames(data)

	} else if (outunit == 'fpkm' || outunit == 'rpkm') {
		# tpm = exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
		# log(tpm) = log(fpkm) - log(sum(fpkm)) + log(1e6)
		#                        ????????????? use nreads, but have to be carefull!
		# log(fpkm) = log(tpm) - log(1e6) + log(nreads)

		samples = colnames(data)
		totalnr = nreadsFromFile(samples, nreads)
		exp = sapply(samples, function(s) {
			exp( log(data[, s]) - log(1e6) + log(totalnr[s, ]) )
		})
		rownames(exp) = rownames(data)

	} else if (outunit == 'cpm') {
		# log(cpm) = log(counts) - log(N) + log(1e6)
		# tpm:
		# rate <- log(counts) - log(effLen)
		# denom <- log(sum(exp(rate)))
		# log(tpm) = rate - denom + log(1e6)
		#
		# log(cpm) - log(tpm) = log(counts) - log(N) - rate + denom
		# log(cpm) - log(tpm) = log(counts) - log(N) - log(counts) + log(effLen) + log(sum(exp(rate)))
		# log(cpm) - log(tpm) = log(effLen) + log(sum(exp(rate))) - log(N)
		#
		# log(sum(exp(rate))) = log(sum(exp(log(counts) - log(effLen))))
		#                     = log(sum(counts/effLen))
		# approximately       = log(sum(counts)/sum(effLen) * length(effLen))
		#                     = log(N) - log(sum(effLen)) + log(length(effLen))
		# 
		# so:
		# log(cpm) - log(tpm) = log(effLen) + log(N) - log(sum(effLen)) + log(length(effLen)) - log(N)
		# log(cpm) = log(tpm) + log(effLen) - log(sum(effLen)) + log(length(effLen))

		glen     = glenFromExon(refexon)
		refgenes = rownames(glen)
		mygenes  = rownames(data)
		outgenes = intersect(refgenes, mygenes)
		if (length(outgenes) < length(mygenes))
			logger('Genes not found in refexon: ', paste(setdiff(mygenes, outgenes), collapse = '\n'), level = 'WARNING')
		data    = data[outgenes, , drop = F]
		glen    = glen[outgenes, , drop = T]

		samples = colnames(data)
		fld     = meanflFromFile(samples, meanfl)
		ngenes  = length(outgenes)

		exp = sapply(samples, function(s) {
			exp( log(data[, s]) + log(glen - fld[s, ] + 1) - log(sum(glen - fld[s, ] + 1)) + log(ngenes) )
		})
		rownames(exp) = outgenes

	} else {
		stop(paste("Don't know yet how to convert fpkm/rpkm to", outunit))
	}

} else if (inunit == 'cpm') {
	
	if (outunit == 'count' || outunit == 'rawcount' || outunit == 'rawcounts' || outunit == 'counts') {
		# log(cpm) = log(counts) - log(sum(counts)) + log(1e6)
		# log(counts) = log(cpm) + log(nreads) - log(1e6)
		# counts = exp( log(cpm) + log(nreads) - log(1e6) )
		samples = colnames(data)
		totalnr = nreadsFromFile(samples, nreads)

		exp = sapply(samples, function(s) {
			exp( log(data[, s]) + log(totalnr[s, ]) - log(1e6) )
		})
		rownames(exp) = rownames(data)
	} else if (outunit == 'fpkm' || outunit == 'rpkm') {
		# log(cpm) = log(counts) - log(sum(counts)) + log(1e6)
		# log(fpkm) = log(counts) + log(1e9) - log(effLen) - log(sum(counts)) 
		# log(fpkm) - log(cpm) = log(1e9) - log(effLen) - log(1e6) = log(1e3) - log(effLen)
		# fpkm = exp( log(cpm) + log(1e3) - log(effLen) )

		glen     = glenFromExon(refexon)
		refgenes = rownames(glen)
		mygenes  = rownames(data)
		outgenes = intersect(refgenes, mygenes)
		if (length(outgenes) < length(mygenes))
			logger('Genes not found in refexon: ', paste(setdiff(mygenes, outgenes), collapse = '\n'), level = 'WARNING')
		data    = data[outgenes, , drop = F]
		glen    = glen[outgenes, , drop = T]

		samples = colnames(data)

		exp = sapply(samples, function(s) {
			exp( log(data[, s]) + log(1e3) - log(glen - fld[s, ] + 1) )
		})
		rownames(exp) = outgenes
	} else if (outunit == 'tpm') {
		# log(cpm) = log(counts) - log(sum(counts)) + log(1e6) [1]
		# tpm:
		# rate <- log(counts) - log(effLen)
		# denom <- log(sum(exp(rate)))
		# log(tpm) = rate - denom + log(1e6)
		# So:
		# log(tpm) = log(counts) - log(effLen) - log(sum(exp(log(counts) - log(effLen)))) + log(1e6)
		# log(tpm) = log(counts) - log(effLen) - log(sum(counts/effLen)) + log(1e6) [2]
		# [2] - [1]:
		# log(tpm) - log(cpm) = log(sum(counts)) - log(effLen) - log(sum(counts/effLen))
		# tpm = exp( log(cpm) + log(sum(counts)) - log(effLen) - log(sum(counts/effLen)) )
		# tpm = exp( log(cpm) + log(nreads)      - log(effLen) - log(?????????????????) )
		# ??? estimated by sum(counts)/sum(effLen)
		# tpm = exp( log(cpm) + log(nreads) - log(effLen) - log(sum(counts)/sum(effLen) * length(effLen)) )
		# tpm = exp( log(cpm) + log(nreads) - log(effLen) - log(sum(counts)) - log(sum(effLen)) + log(length(effLen)) )
		# tpm = exp( log(cpm) - log(effLen) - log(sum(effLen)) + log(length(effLen)) )
		glen     = glenFromExon(refexon)
		refgenes = rownames(glen)
		mygenes  = rownames(data)
		outgenes = intersect(refgenes, mygenes)
		if (length(outgenes) < length(mygenes))
			logger('Genes not found in refexon: ', paste(setdiff(mygenes, outgenes), collapse = '\n'), level = 'WARNING')
		data    = data[outgenes, , drop = F]
		glen    = glen[outgenes, , drop = T]

		samples = colnames(data)
		ngenes  = length(outgenes)

		exp = sapply(samples, function(s) {
			exp( log(data[, s]) - log(glen - fld[s, ] + 1) - log(sum(glen - fld[s, ] + 1)) + log(ngenes) )
		})
		rownames(exp) = outgenes
	}

}

if (!is.null(outform)) exp = outform(exp)

write.table (round(exp, 3), outfile, quote=F, row.names=T, col.names=T, sep="\t")

