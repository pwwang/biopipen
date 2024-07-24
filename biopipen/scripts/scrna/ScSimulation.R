{{ biopipen_dir | joinpaths: "utils", "misc.R" | source_r }}

library(rlang)
library(splatter)
library(scater)

# Load template variables
seed <- {{ in.seed | r }}
outfile <- {{ out.outfile | r }}
ngenes <- {{ envs.ngenes | r }}
ncells <- {{ envs.ncells | r }}
nspikes <- {{ envs.nspikes | r }}
outtype <- {{ envs.outtype | r }}
method <- {{ envs.method | r }}
user_params <- {{ envs.params | r: todot="-" }}

log_info("Generating simulation parameters ...")

seed <- seed %||% 1
if (length(seed) > 1) {
    log_warn("- multiple seeds provided, using the first one")
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


log_info("Saving simulation parameters to file ...")

sim <- splatSimulate(params, method = method, verbose = TRUE)

outtype <- tolower(outtype)
if (outtype == "sce") outtype <- "singlecellexperiment"

if (outtype == "singlecellexperiment") {
    log_info("Saving simulation to file ...")
    saveRDS(sim, file = outfile)
} else {
    log_info("Converting simulation to Seurat object ...")
    cnts <- SingleCellExperiment::counts(sim)
    sobj <- Seurat::CreateSeuratObject(counts = cnts, project = proj)
    rm(sim)
    rm(cnts)
    gc()

    log_info("Saving simulation to file ...")
    saveRDS(sobj, file = outfile)
}
