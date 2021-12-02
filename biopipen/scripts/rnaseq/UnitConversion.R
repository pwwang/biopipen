infile = {{in.infile | r}}
outfile = {{out.outfile | r}}
infmt = {{envs.infmt | r}}
inunit = {{envs.inunit | r}}
outunit = {{envs.outunit | r}}
refexon = {{envs.refexon | r}}
inlog2p = {{envs.inlog2p | r}}
outlog2p = {{envs.outlog2p | r}}

if (infmt == "rds") {
    indata = readRDS(infile)
} else if (endsWith(infile, ".gz")) {
    indata = read.table(infile, header=T, row.names=NULL, sep="\t", check.names = F)
    genes = make.unique(indata[, 1])
    indata = indata[, -1]
    rownames(indata) = genes
}


glenFromExon = function(exonfile, x) {
	gff  = read.table(exonfile, header = F, row.names = NULL)
	# V4: start, V5: end, V10: gene name
	glen = aggregate(V5-V4+1 ~ V10, gff, sum)
	genes = glen[,1]
	glen = glen[,-1,drop=F]
	rownames(glen) = genes

	mygenes  = rownames(x)
	outgenes = intersect(genes, mygenes)
	if (length(outgenes) < length(mygenes))
		warning('Genes not found in refexon: ', paste(setdiff(mygenes, outgenes)))

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

count2tpm = function(x) {
	glen = glenFromExon(refexon, x)
	x = x[rownames(glen), , drop = F]
	fld  = meanflFromFile(samples, meanfl)

	# see: https://gist.github.com/slowkow/c6ab0348747f86e2748b
	expr = as.data.frame(sapply(samples, function(s){
		rate  = log(x[, s]) - log(glen - fld[s, ] + 1)
		denom = log(sum(exp(rate)))
		exp(rate - denom + log(1e6))
	}))
	colnames(expr) = colnames(x)
	rownames(expr) = rownames(x)
	expr
}

if (inunit %in% c('count', 'counts', 'rawcount', 'rawcounts')) {
    inunit = "count"
}


convert = function(data, inunit, outunit) {
    func = get(paste0(inunit, "2", outunit))
    func(data)
}

convert(indata, inunit, outunit)
