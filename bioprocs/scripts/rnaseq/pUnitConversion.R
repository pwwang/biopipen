{{rimport}}('__init__.r')

infile  = {{i.infile | R}}
outfile = {{o.outfile | R}}
inunit  = {{args.inunit | R}}
outunit = {{args.outunit | R}}
meanfl  = {{args.meanfl | R}}
nreads  = {{args.nreads | R}}
refexon = {{args.refexon | R}}
inform  = {{args.inform | lambda x: 'NULL' if not x else x}}
outform = {{args.outform | lambda x: 'NULL' if not x else x}}

data    = read.table.inopts (infile, list(cnames = TRUE), dup = 'drop')
samples = colnames(data)
if (!is.null(inform)) {
	data = inform(data)
}

glenFromExon = function(exonfile, data) {
	gff  = read.table(exonfile, header = F, row.names = NULL)
	# V4: start, V5: end, V10: gene name
	glen = aggregate(V5-V4+1 ~ V10, gff, sum)
	genes = glen[,1]
	glen = glen[,-1,drop=F]
	rownames(glen) = genes

	mygenes  = rownames(data)
	outgenes = intersect(refgenes, mygenes)
	if (length(outgenes) < length(mygenes))
		logger('Genes not found in refexon: ', paste(setdiff(mygenes, outgenes), collapse = ','), level = 'WARNING')

	glen[outgenes, , drop = FALSE]
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

count2cpm = function(data) {
	edgeR::cpm(data)
}

count2fpkm = function(data) {
	# may lose some genes
	glen = glenFromExon(refexon)
	data = data[rownames(glen), , drop = F]
	dge  = edgeR::DGEList(counts=data)

	dge$genes$Length = glen
	edgeR::rpkm(dge)
}

count2fpkmuq = function(data) {
	# may lose some genes
	glen = glenFromExon(refexon)
	data = data[rownames(glen), , drop = FALSE]

	fld  = meanflFromFile(samples, meanfl)
	expr = sapply(samples, function(s){
		RC75 = quantile(data[, s], .75)
		exp( log(data[, s]) + log(1e9) - log(glen - fld[s, ] + 1) - log(RC75) )
	})
	rownames(expr) = rownames(data)
	expr
}

count2tpm = function(data) {
	glen = glenFromExon(refexon)
	data = data[rownames(glen), , drop = F]
	fld  = meanflFromFile(samples, meanfl)

	# see: https://gist.github.com/slowkow/c6ab0348747f86e2748b
	expr     = sapply(samples, function(s){
		rate  = log(data[, s]) - log(glen - fld[s, ] + 1)
		denom = log(sum(exp(rate)))
		exp(rate - denom + log(1e6))
	})
	rownames(expr) = rownames(data)
	expr
}

count2tmm = function(data) {
	dge = edgeR::DGEList(counts=data)
	dge = edgeR::calcNormFactors(dge, method = "TMM")
	edgeR::cpm(dge)
}

fpkm2count = function(data) {
	glen    = glenFromExon(refexon)
	data    = data[rownames(glen), , drop = F]
	fld     = meanflFromFile(samples, meanfl)
	totalnr = nreadsFromFile(samples, nreads)

	expr    = sapply(samples, function(s){
		N = totalnr[s, ]
		exp( log(data[, s]) + log(N) + log(glen - fld[s, ] + 1) - log(1e9) )
	})
	rownames(expr) = rownames(data)
	expr
}

fpkm2tpm = function(data) {
	expr = sapply(samples, function(s) {
		exp( log(data[, s]) - log(sum(data[, s])) + log(1e6) )
	})
	rownames(expr) = rownames(data)
	expr
}

fpkm2cpm = function(data) {
	glen = glenFromExon(refexon)
	data = data[rownames(glen), , drop = F]
	expr = sapply(samples, function(s) {
		exp( log(data[, s]) - log(1e3) - log(glen - fld[s, ] + 1) )
	})
	rownames(expr) = rownames(data)
	expr
}

tpm2count = function(data) {
	totalnr = nreadsFromFile(samples, nreads)
	ngenes  = nrow(data)

	expr     = sapply(samples, function(s){
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
	rownames(expr) = rownames(data)
	expr
}

tpm2fpkm = function(data) {
	totalnr = nreadsFromFile(samples, nreads)
	expr = sapply(samples, function(s) {
		exp( log(data[, s]) - log(1e6) + log(totalnr[s, ]) )
	})
	rownames(expr) = rownames(data)
	expr
}

tpm2cpm = function(data) {
	glen   = glenFromExon(refexon)
	data   = data[rownames(glen), , drop = F]
	fld    = meanflFromFile(samples, meanfl)
	ngenes = length(outgenes)

	expr = sapply(samples, function(s) {
		exp( log(data[, s]) + log(glen - fld[s, ] + 1) - log(sum(glen - fld[s, ] + 1)) + log(ngenes) )
	})
	rownames(expr) = rownames(data)
	expr
}

cpm2count = function(data) {
	totalnr = nreadsFromFile(samples, nreads)

	expr = sapply(samples, function(s) {
		exp( log(data[, s]) + log(totalnr[s, ]) - log(1e6) )
	})
	rownames(expr) = rownames(data)
	expr
}

cpm2fpkm = function(data) {
	glen = glenFromExon(refexon)
	data = data[rownames(glen), , drop = F]
	expr = sapply(samples, function(s) {
		exp( log(data[, s]) + log(1e3) - log(glen - fld[s, ] + 1) )
	})
	rownames(expr) = rownames(data)
	expr
}

cpm2tpm = function(data) {
	glen   = glenFromExon(refexon)
	data   = data[rownames(glen), , drop = F]
	ngenes = nrow(glen)
	expr   = sapply(samples, function(s) {
		exp( log(data[, s]) - log(glen - fld[s, ] + 1) - log(sum(glen - fld[s, ] + 1)) + log(ngenes) )
	})
	rownames(expr) = rownames(data)
	expr
}

is.count  = function(unit) {unit %in% c('count', 'counts', 'rawcount', 'rawcounts')}
is.cpm    = function(unit) {unit == 'cpm'}
is.fpkm   = function(unit) {unit %in% c('fpkm', 'rpkm')}
is.fpkmuq = function(unit) {unit %in% c('fpkmuq', 'rpkmuq', 'fpkm-uq', 'rpkm-uq')}
is.tpm    = function(unit) {unit == 'tpm'}
is.tmm    = function(unit) {unit == 'tmm'}


convert = function(data, inunit, outunit) {
	if (is.count(inunit) && is.cpm(outunit))
		return (count2cpm(data))
	if (is.count(inunit) && is.fpkm(outunit))
		return (count2fpkm(data))
	if (is.count(inunit) && is.fpkmuq(outunit))
		return (count2fpkmuq(data))
	if (is.count(inunit) && is.tpm(outunit))
		return (count2tpm(data))
	if (is.count(inunit) && is.tmm(outunit))
		return (count2tmm(data))
	if (is.fpkm(inunit) && is.count(outunit))
		return (fpkm2count(data))
	if (is.fpkm(inunit) && is.tpm(outunit))
		return (fpkm2tpm(data))
	if (is.fpkm(inunit) && is.cpm(outunit))
		return (fpkm2cpm(data))
	if (is.tpm(inunit) && is.count(outunit))
		return (tpm2count(data))
	if (is.tpm(inunit) && is.fpkm(outunit))
		return (tpm2fpkm(data))
	if (is.tpm(inunit) && is.cpm(outunit))
		return (tpm2cpm(data))
	if (is.cpm(inunit) && is.count(outunit))
		return (cpm2count(data))
	if (is.cpm(inunit) && is.fpkm(outunit))
		return (cpm2fpkm(data))
	if (is.cpm(inunit) && is.tpm(outunit))
		return (cpm2tpm(data))
	stop(paste('Unsupported conversion from', inunit, 'to', outunit, '.'))
}

expr = convert(data, inunit, outunit)
if (!is.null(outform)) {
	expr = outform(expr)
}

write.table(round(expr, 3), outfile, quote=FALSE, row.names=TRUE, col.names=TRUE, sep="\t")
