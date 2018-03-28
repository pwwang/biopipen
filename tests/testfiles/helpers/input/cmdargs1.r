library(reticulate)
bioprocs = import('bioprocs')
source(file.path(bioprocs$UTILS, 'helpers.r'))

params = list()
dash   = 'auto'
equal  = 'auto'
cat(cmdargs(params, dash, equal))