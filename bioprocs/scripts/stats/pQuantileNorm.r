library(methods)
library(preprocessCore)
{{rimport}}('__init__.r')

infile  = {{i.infile | R}}
outfile = {{o.outfile | R}}
inopts  = {{args.inopts | R}}

indata = read.table.inopts(infile, inopts)
rnames = ifelse(inopts$rnames, rownames(indata), NULL)
cnames = ifelse(inopts$cnames, colnames(indata), NULL)

indata = normalize.quantiles(as.matrix(indata))
if (!is.null(rnames)) rownames(indata) = rnames
if (!is.null(cnames)) colnames(indata) = cnames

write.table(indata, outfile, col.names = inopts$cnames, row.names = inopts$rnames, sep = "\t", quote = F)