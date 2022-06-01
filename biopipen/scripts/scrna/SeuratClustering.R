library(Seurat)
library(tidyr)
library(dplyr)

set.seed(8525)

srtfile = {{in.srtobj | quote}}
rdsfile = {{out.rdsfile | quote}}
envs = {{envs | r}}

sobj = readRDS(srtfile)

sobj = FindNeighbors(sobj)
envs$FindClusters$object = sobj
sobj = do.call(FindClusters, envs$FindClusters)
saveRDS(sobj, file = rdsfile)
