
library(Seurat)
library(biopipen.utils)

set.seed(8525)

srtfile <- {{in.srtobj | r}}
outfile <- {{out.outfile | r}}
joboutdir <- {{job.outdir | r}}
RunPCAArgs <- {{envs.RunPCA | r: todot="-"}}
FindNeighborsArgs <- {{envs.FindNeighbors | r: todot="-"}}
FindClustersArgs <- {{envs.FindClusters | r: todot="-"}}
RunUMAPArgs <- {{envs.RunUMAP | r: todot="-"}}
cache <- {{envs.cache | r}}
ncores <- {{envs.ncores | r}}

log <- get_logger()

# options(str = strOptions(vec.len = 5, digits.d = 5))
options(future.globals.maxSize = 80000 * 1024^2)
plan(strategy = "multicore", workers = ncores)

log$info("Reading Seurat object ...")
sobj <- read_obj(srtfile)

if (isTRUE(cache)) { cache <- joboutdir }

sobj <- RunSeuratClustering(
    sobj,
    RunPCAArgs = RunPCAArgs,
    RunUMAPArgs = RunUMAPArgs,
    FindNeighborsArgs = FindNeighborsArgs,
    FindClustersArgs = FindClustersArgs,
    log = log,
    cache = cache
)

log$info("Saving results ...")
save_obj(sobj, file = outfile)
