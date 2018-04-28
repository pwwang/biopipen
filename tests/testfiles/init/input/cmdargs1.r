library(reticulate)
bioprocs = import('bioprocs')
source(file.path(bioprocs$UTILS, '__init__.r'))

params = list()
dash   = 'auto'
equal  = 'auto'
cat(cmdargs(params, dash, equal))