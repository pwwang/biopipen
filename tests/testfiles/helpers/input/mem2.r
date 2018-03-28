library(reticulate)
bioprocs = import('bioprocs')
source(file.path(bioprocs$UTILS, 'helpers.r'))

args = commandArgs(trailingOnly=TRUE)
cat(mem2(args[1], args[2]))