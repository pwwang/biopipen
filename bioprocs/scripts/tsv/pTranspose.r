{{rimport}}('__init__.r')

infile  = {{i.infile | R}}
outfile = {{o.outfile | R}}
inopts  = {{args.inopts | R}}

indata = t(read.table.inopts(infile, inopts))
write.table(indata, outfile, quote = F, sep = "\t", row.names = inopts$cnames, col.names = inopts$rnames)