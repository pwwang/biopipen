library(reticulate)
bioprocs = import('bioprocs')
source(file.path(bioprocs$UTILS, '__init__.r'))

params = list(c=T, de='de fg', b=2, a=1)
dash   = '--'
equal  = ' '
cat(cmdargs(params, dash, equal))