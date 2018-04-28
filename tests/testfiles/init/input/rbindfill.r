library(reticulate)
bioprocs = import('bioprocs')
source(file.path(bioprocs$UTILS, '__init__.r'))

args = commandArgs(trailingOnly=TRUE)
infile1 = args[1]
infile2 = args[2]
outfile = args[3]

data1 = read.table(infile1, header = T, row.names = 1, check.names = F)
data2 = read.table(infile2, header = T, row.names = 1, check.names = F)

d = rbindfill(data1, data2)
d = d[order(rownames(d)), order(colnames(d))]
write.table(d, outfile, sep = '\t', quote = F, row.names = T, col.names = T)