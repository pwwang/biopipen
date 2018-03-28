library(reticulate)
bioprocs = import('bioprocs')
source(file.path(bioprocs$UTILS, 'sampleinfo.r'))

args = commandArgs(trailingOnly=TRUE)
sfile = args[1]
nrow  = args[2]
ncol  = args[3]
si = SampleInfo(sfile)
d  = si$dataframe(si$data)

stopifnot(nrow(d) == nrow)
stopifnot(ncol(d) == ncol)