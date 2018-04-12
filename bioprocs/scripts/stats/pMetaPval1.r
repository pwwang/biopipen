library("methods")
library("metap")

infile  = {{in.infile | quote}}
outfile = {{out.outfile | quote}}

inopts.cnames = {{args.inopts.cnames | R}}
inopts.pcol   = {{args.inopts.pcol | R}}
outopts.head  = {{args.outopts.head | R}}
outopts.ponly = {{args.outopts.ponly | R}}
method = {{args.method | quote}}

indata = read.table(infile, header = inopts.cnames, row.names = NULL, check.names = F)
rnames = unique(as.vector(indata[, 1]))
if (inopts.pcol < 0) inopts.pcol = ncol(indata) + inopts.pcol + 1

ret = NULL
for (rname in rnames) {
	pvals               = as.vector(indata[which(indata[,1] == rname), inopts.pcol])
	pvals[is.na(pvals)] = 1
	pvals[pvals > 1]    = 1
	pvals[pvals <= 0]   = 1e-100
	if(length(pvals) == 1) {
		mp1 = if (method == 'logitp')
			list(t = 0, df = 0, p = pvals[1], validp = pvals[1]) else if (method == 'sumz')
			list(z = 0, p = pvals[1], validp = pvals[1]) else if (method == 'votep')
			list(p = pvals[1], pos = 0, neg = 0, alpha = 0, validp = pvals[1]) else if (method == 'sump')
			list(p = pvals[1], conservativep = pvals[1], validp = pvals[1]) else if (method == 'meanp')
			list(z = 0, p = pvals[1], validp = pvals[1]) else if (method == 'wilkinsonp')
			list(p = pvals[1], pr = 0, r = 0, critp = pvals[1], alpha = 0, validp = pvals[1])
	} else {
		mp1 = {{args.method}}(pvals)
	}
	mp2    = unlist(mp1)
	mnames = names(mp1)
	mnames = if ("weights" %in% mnames) mnames[1:length(mnames)-2] else mnames[1:length(mnames)-1]
	mmp    = matrix(mp2[1:length(mnames)], nrow = 1)
	rownames(mmp) = rname
	colnames(mmp) = mnames
	mmp    = cbind(mmp, validp = paste(mp1$validp, collapse=';'))
	ret    = if (is.null(ret)) mmp else rbind(ret, mmp)
}
ret = if (outopts.ponly) ret[order(as.numeric(ret[, 'p'])), 'p', drop=F] else ret[order(as.numeric(ret[, 'p'])),,drop=F]
write.table(ret, outfile, sep="\t", quote=F, col.names=outopts.head)
