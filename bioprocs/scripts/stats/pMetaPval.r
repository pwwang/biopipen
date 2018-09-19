library("methods")
library("metap")
{{rimport}}('__init__.r')

infiles = Sys.glob(file.path({{i.indir | quote}}, {{args.pattern | quote}}))
lenfile = length(infiles)

outfile = {{o.outfile | quote}}

inopts.cnames = {{args.inopts.cnames | lambda x: x if isinstance(x, list) else [x] | Rvec}}
inopts.pcol   = {{args.inopts.pcol   | lambda x: x if isinstance(x, list) else [x] | Rvec}}
inopts.cnames = as.logical(inopts.cnames)

outopts.head  = {{args.outopts.head | R}}
outopts.ponly = {{args.outopts.ponly | R}}

if (length(inopts.cnames) == 1) {
	cnames = rep(inopts.cnames, lenfile)
} else {
	cnames = rep(T, lenfile)
	for (i in 1:length(inopts.cnames)) {
		cnames[i] = inopts.cnames[i]
	}
}

if (length(inopts.pcol) == 1) {
	pcol = rep(inopts.pcol, lenfile)
} else {
	cnames = rep(-1, lenfile)
	for (i in 1:length(inopts.pcol)) {
		cnames[i] = inopts.pcol[i]
	}
}

pvals = NULL
for (i in 1:lenfile) {
	infile = infiles[i]
	indata = read.table(infile, sep="\t", header=cnames[i], row.names=NULL, check.names=F)
	rnames = make.unique(as.vector(indata[,1]))
	indata[,1] = NULL
	rownames(indata) = rnames
	col    = pcol[i]
	if (col<0) col = ncol(indata) + col + 1
	indata = indata[, col, drop=F]
	colnames(indata) = paste(colnames(indata), tools::file_path_sans_ext(basename(infile)), sep='_')
	if (is.null(pvals)) {
		pvals = indata
	} else {
		pvals = cbindfill(pvals, indata)
	}
}
pvals[is.na(pvals)] = 1
pvals[pvals > 1]    = 1
pvals[pvals <= 0]   = 1e-100

ret     = apply(pvals, 1, {{args.method}})
ret2    = matrix(unlist(ret), byrow=T, nrow=nrow(pvals))
cnames  = colnames(pvals)
rnames  = rownames(pvals)
rnames2 = names(ret[[rnames[1]]])
rnames2 = if ("weights" %in% rnames) rnames2[1:length(rnames2)-2] else rnames2[1:length(rnames2)-1]
rownames(ret2) = rnames
colnames(ret2) = c(rnames2, cnames)

if (outopts.ponly) {
	ret2 = ret2[order(ret2[, 'p']), 'p', drop=F]
} else {
	ret2 = ret2[order(ret2[, 'p']),,drop=F]
}
write.table(ret2, outfile, sep="\t", quote=F, col.names=outopts.head)
