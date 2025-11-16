library(Seurat)
library(rlang)
library(biopipen.utils)

set.seed(8525)

sobjfile = {{in.sobjfile | r}}
outfile = {{out.outfile | r}}
joboutdir = {{job.outdir | r}}
use = {{envs.use | r}}
ident = {{envs.ident | r}}
ref = {{envs.ref | r}}
refnorm = {{envs.refnorm | r}}
ncores = {{envs.ncores | r}}
split_by = {{envs.split_by | r}}
mutaters = {{envs.mutaters | r}}
skip_if_normalized = {{envs.skip_if_normalized | r}}
sctransform_args = {{envs.SCTransform | r: todot="-"}}
normalizedata_args = {{envs.NormalizeData | r: todot="-"}}
findtransferanchors_args = {{envs.FindTransferAnchors | r: todot="-"}}
mapquery_args = {{envs.MapQuery | r: todot="-"}}
cache = {{envs.cache | r}}
plots = {{envs.plots | r}}

log <- get_logger()
reporter <- get_reporter()

options(future.globals.maxSize = Inf)
options(future.rng.onMisuse="ignore")
options(Seurat.object.assay.version = "v5")

# See if we have a reference
if (is.null(ref)) {
    stop("No reference provided (envs.ref)")
}

if (is.null(use)) {
    stop("No use provided (envs.use), don't know which column to transfer as cluster")
}

outdir = dirname(outfile)
if (isTRUE(cache)) {
    cache = joboutdir
}
if (is.null(split_by)) {
    future::plan(strategy = "multicore", workers = ncores)
}

log$info("Loading reference ...")
if (endsWith(ref, ".rds") || endsWith(ref, ".RDS") || endsWith(ref, ".qs") || endsWith(ref, ".qs2")) {
    reference <- read_obj(ref)
} else if (endsWith(ref, ".h5seurat") || endsWith(ref, ".H5Seurat")) {
    reference <- SeuratDisk::LoadH5Seurat(ref)
} else {
    stop("Reference file must be .qs, .qs2, .rds, .RDS, .h5seurat or .H5Seurat")
}
reference <- tryCatch(JoinLayers(reference), error = function(e) {reference})
Idents(reference) <- reference@meta.data[[use]]

log$info("Loading query data ...")
sobj <- read_obj(sobjfile)

sobj <- RunSeuratMap2Ref(
    object = sobj, ref = reference, use = use,
    ident = ident, refnorm = refnorm, skip_if_normalized = skip_if_normalized,
    split_by = split_by, ncores = ncores,
    SCTransformArgs = sctransform_args,
    NormalizeDataArgs = normalizedata_args,
    FindTransferAnchorsArgs = findtransferanchors_args,
    MapQueryArgs = mapquery_args,
    log = log, cache = cache
)

# Save
gc()
log$info("Saving result ...")
save_obj(sobj, file = outfile)


### Plotting
log$info("Plotting features ...")
for (name in names(plots)) {
    if (is.null(plots[[name]])) {
        next
    }
    log$info("- {name} ...")
    plots[[name]]$features <- gsub("{use}", use, plots[[name]]$features, fixed = TRUE)
    plots[[name]]$features <- gsub("{ident}", ident, plots[[name]]$features, fixed = TRUE)

    plots[[name]]$devpars <- plots[[name]]$devpars %||% list()
    plots[[name]]$devpars$res <- plots[[name]]$devpars$res %||% 100
    plots[[name]]$devpars$width <- plots[[name]]$devpars$width %||% 1200
    plots[[name]]$devpars$height <- plots[[name]]$devpars$height %||% 720
    plots[[name]]$more_formats <- plots[[name]]$more_formats %||% character()
    plots[[name]]$save_code <- FALSE
    plots[[name]]$descr <- plots[[name]]$descr %||% name
    extract_vars(plots[[name]], "devpars", "more_formats", "save_code", "descr")

    plot_fn <- gglogger::register(VizSeuratMap2Ref)
    p <- do_call(plot_fn, c(list(query = sobj, ref = reference), plots[[name]]))
    prefix <- file.path(outdir, paste0(slugify(name), ".map2ref"))
    save_plot(p, prefix, devpars, formats = c("png", more_formats))

    reporter$add(
        reporter$image(prefix, more_formats, save_code = FALSE, kind = "image"),
        h1 = name
    )
}

reporter$save(joboutdir)
