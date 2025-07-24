library(rlang)
library(splatter)
library(scater)
library(biopipen.utils)

# Load template variables
seed <- {{ in.seed | r }}
outfile <- {{ out.outfile | r }}
ngenes <- {{ envs.ngenes | r }}
ncells <- {{ envs.ncells | r }}
nspikes <- {{ envs.nspikes | r }}
outtype <- {{ envs.outtype | r }}
method <- {{ envs.method | r }}
user_params <- {{ envs.params | r: todot="-" }}

log <- get_logger()

log$info("Generating simulation parameters ...")

seed <- seed %||% 1
if (length(seed) > 1) {
    log$warn("- multiple seeds provided, using the first one")
    seed <- seed[1]
}
if (is.character(seed)) {
    library(digest)
    proj <- seed
    seed <- digest2int(seed)
} else {
    proj <- paste0("S", seed)
}

set.seed(seed)
mock_sce_params <- list()
if (!is.null(ngenes)) mock_sce_params$ngenes <- ngenes
if (!is.null(ncells)) mock_sce_params$ncells <- ncells
if (!is.null(nspikes)) mock_sce_params$nspikes <- nspikes
sce <- do.call(mockSCE, mock_sce_params)
params <- splatEstimate(sce)
user_params$seed <- seed
user_params$object = params
do_call(setParams, user_params)


log$info("Saving simulation parameters to file ...")

sim <- splatSimulate(params, method = method, verbose = TRUE)

outtype <- tolower(outtype)
if (outtype == "sce") outtype <- "singlecellexperiment"

if (outtype == "singlecellexperiment") {
    log$info("Saving simulation to file ...")
    save_obj(sim, file = outfile)
} else {
    log$info("Converting simulation to Seurat object ...")
    cnts <- SingleCellExperiment::counts(sim)
    sobj <- Seurat::CreateSeuratObject(counts = cnts, project = proj)
    rm(sim)
    rm(cnts)
    gc()

    log$info("Saving simulation to file ...")
    save_obj(sobj, file = outfile)
}
