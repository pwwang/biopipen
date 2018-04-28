library(reticulate)
bioprocs = import('bioprocs')
source(file.path(bioprocs$UTILS, '__init__.r'))

runcmd('cmdnotexists')