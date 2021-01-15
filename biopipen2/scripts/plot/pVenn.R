{{'plot.r', '__init__.r' | rimport}}

infile   = {{i.infile | R}}
metafile = {{i.metafile | R}}
outfile  = {{o.outfile | R}}
inopts   = {{args.inopts | R}}
intype   = {{args.intype | R}}

inopts$fill = TRUE
data = read.table.inopts(infile, inopts)
if (intype == 'raw') {
	allelements = unique(as.vector(as.matrix(data)))
	allelements = allelements[allelements!=""]
	data = as.data.frame(apply(data, 2, function(x) as.integer(allelements %in% x)))
}

tool    = {{args.tool | R}}
ncols   = ncol(data)
cnames  = colnames(data)
devpars = {{args.devpars | R}}
# use real venn plot
if ((ncols <= 3 && tool == 'auto') || tool == 'venn') {
	params = {{args.params | R}}
	plot.venn(data, outfile, params, devpars)
} else {
	metadata = data.frame(counts = colSums(data), ratio = colSums(data)/nrow(data))
	if (metafile != "") {
		exmeta   = read.table(metafile, sep = "\t", header = TRUE, row.names = 1, check.names = F)
		metadata = cbind.fill(metadata, exmeta)
	}
	metadata = cbind(metadata, sets = rownames(metadata))
	missingcols = rownames(metadata[metadata$counts == 0,,drop=F])
	if (length(missingcols) > 0) {
		log2pyppl('No data in columns:', missingcols)
	}
	metadata = metadata[metadata$counts>0, , drop = F]
	data = data[, rownames(metadata), drop = F]
	params = {{args.params | R}}
	plot.upset(data, outfile, params, devpars)
}
