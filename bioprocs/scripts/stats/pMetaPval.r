library("metap")
{{cbindfill}}

infiles = list.files(path={{in.indir | quote}}, pattern = {{args.pattern | quote}}, full.names=T)
lenfile = length(infiles)
header  = {{args.header | R}}
pcol    = {{args.pcol | R}}
if (length(header) < lenfile) {
	header = c(header, rep(T, lenfile - length(header)))
}
if (length(pcol) < lenfile) {
	pcol = c(pcol, rep(-1, lenfile - length(pcol)))
}

pvals = NULL
for (i in 1:lenfile) {
	infile = infiles[i]
	indata = read.table(infile, sep="\t", header=header[i], row.names=NULL, check.names=F)
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

ret2 = {% if args.poutonly %}ret2[order(ret2[, 'p']), 'p', drop=F]{% else %}ret2[order(ret2[, 'p']),,drop=F]{% endif %}
write.table(ret2, {{out.outfile | quote}}, sep="\t", quote=F, col.names={{args.outheader | R}})
