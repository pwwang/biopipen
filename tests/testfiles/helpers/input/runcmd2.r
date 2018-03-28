library(reticulate)
bioprocs = import('bioprocs')
source(file.path(bioprocs$UTILS, 'helpers.r'))

runcmd('cmdnotexists')