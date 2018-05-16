{{rimport}}('plot.r', '__init__.r')

rnames  = {{args.rnames | R}}
infile  = {{in.infile | R}}
outfile = {{out.outfile | R}}

data = read.table.nodup ({{in.infile | quote}}, sep="\t", header = TRUE, row.names = rnames, check.names = F)

tool    = {{args.tool | R}}
ncols   = ncol(data)
cnames  = colnames(data)
params  = {{args.params | R}}
devpars = {{args.devpars | R}}
# use real venn plot
if ((ncols <= 3 && tool == 'auto') || tool == 'venn') {
	plot.venn(data, outfile, params, devpars)
} else {
	plot.upset(data, outfile, params, devpars)
}
