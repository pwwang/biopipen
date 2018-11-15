library(philentropy)

infile  = {{i.infile | R}}
outfile = {{o.outfile | R}}
inopts  = {{args.inopts | R}}
method  = {{args.method | R}}
sim     = {{args.sim | R}}
na      = {{args.na | R}}

indata = read.table(infile, header = inopts$cnames, row.names = if(inopts$rnames) 1 else NULL, sep = "\t", check.names = F)
indata[is.na(indata)] = na

d = distance(indata, method = method, test.na = F)
if (sim) d = 1 - d
rownames(d) = rownames(indata)
colnames(d) = rownames(indata)
write.table(round(d, 3), outfile, col.names = T, row.names = T, sep = "\t", quote = F)
