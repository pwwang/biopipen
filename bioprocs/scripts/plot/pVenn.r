{{rimport}}('plot.r', '__init__.r')

rnames   = {{args.rnames | R}}
infile   = {{in.infile | R}}
metafile = {{in.metafile | R}}
outfile  = {{out.outfile | R}}

data = read.table.nodup (infile, sep="\t", header = TRUE, row.names = rnames, check.names = F)

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
	params = {{args.params | R}}
	plot.upset(data, outfile, params, devpars)
}
