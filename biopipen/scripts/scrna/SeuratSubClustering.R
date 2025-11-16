
library(Seurat)
library(future)
library(biopipen.utils)

set.seed(8525)

srtfile <- {{in.srtobj | r}}
outfile <- {{out.outfile | r}}
ncores <- {{envs.ncores | r}}
mutaters <- {{envs.mutaters | r}}
subset <- {{envs.subset | r}}
cache <- {{envs.cache | r}}
RunPCAArgs <- {{envs.RunPCA | r: todot = "-"}}
RunUMAPArgs <- {{envs.RunUMAP | r: todot = "-"}}
FindNeighborsArgs <- {{envs.FindNeighbors | r: todot = "-"}}
FindClustersArgs <- {{envs.FindClusters | r: todot = "-"}}
cases <- {{envs.cases | r}}

options(future.globals.maxSize = Inf)
plan(strategy = "multicore", workers = ncores)

log <- get_logger()

cases <- expand_cases(cases, defaults = list(
    RunPCA = RunPCAArgs,
    RunUMAP = RunUMAPArgs,
    FindNeighbors = FindNeighborsArgs,
    FindClusters = FindClustersArgs,
    subset = subset
))

if (isTRUE(cache)) {}

log$info("Reading Seurat object ...")
object <- read_obj(srtfile)

if (!is.null(mutaters) && length(mutaters) > 0) {
    log$info("Mutating meta data ...")
    object@meta.data <- mutate(
        object@meta.data,
        !!!lapply(mutaters, parse_expr)
    )
}

for (name in names(cases)) {
    case <- cases[[name]]
    log$info("Processing case '{name}' ...")

    object <- RunSeuratSubClustering(
        object = object,
        subset = case$subset,
        name = name,
        RunPCAArgs = case$RunPCAArgs,
        RunUMAPArgs = case$RunUMAPArgs,
        FindNeighborsArgs = case$FindNeighborsArgs,
        FindClustersArgs = case$FindClustersArgs,
        log = log,
        cache = cache
    )
}

log$info("Saving results ...")
biopipen.utils::save_obj(object, file = outfile)
