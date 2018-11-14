library(philentropy)
infile  = {{i.infile | R}}
outfile = {{o.outfile | R}}
inopts  = {{args.inopts | R}}
method  = {{args.method | R}}
sim     = {{args.sim | R}}

indata = read.table(infile, header = inopts$cnames, row.names = if(inopts$cnames) 1 else NULL, sep = "\t", check.names = F)

d = as.matrix(dist(indata, method = method))
if (sim) d = 1 - d
write.table(d, outfile, col.names = T, row.names = T, sep = "\t", quote = F)
