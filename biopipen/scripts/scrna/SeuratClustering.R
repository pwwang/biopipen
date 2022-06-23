source("{{biopipen_dir}}/utils/misc.R")

library(Seurat)
library(tidyr)
library(dplyr)

set.seed(8525)

srtfile = {{in.srtobj | quote}}
rdsfile = {{out.rdsfile | quote}}
envs = {{envs | r}}

.expand_dims = function(args, name = "dims") {
    # Expand dims from 30 to 1:30
    if (is.numeric(args[[name]]) && length(args[[name]] == 1)) {
        args[[name]] = 1:args[[name]]
    }
    args
}
envs$FindNeighbors = .expand_dims(envs$FindNeighbors)

sobj = readRDS(srtfile)

envs$FindNeighbors$object = sobj
sobj = do_call(FindNeighbors, envs$FindNeighbors)

envs$FindClusters$object = sobj
sobj = do_call(FindClusters, envs$FindClusters)
saveRDS(sobj, file = rdsfile)
