library(rlang)
library(dplyr)
library(Seurat)
library(biopipen.utils)

sobjfile <- {{in.srtobj | r}}
outfile <- {{out.rdsfile | r}}
defaults <- {{envs.defaults | r}}
ncores <- {{envs.ncores | r}}
modules <- {{envs.modules | r}}
post_mutaters <- {{envs.post_mutaters | r}}

qs2::qopt("nthreads", value = ncores)
log <- get_logger()

# load seurat object
log$info("Loading Seurat object ...")
sobj <- read_obj(sobjfile)

# calculate module scores
sobj <- do_call(RunModuleScoring, c(list(object = sobj, modules = modules), defaults))

if (!is.null(post_mutaters) && length(post_mutaters) > 0) {
    log$info("Applying post mutaters ...")
    # MutateSeuratMeta() applies the same dplyr::mutate() but also understands
    # the `<column>:ident` key, so a "<column>" can become the object identity.
    # It has to happen here rather than in CellTypeAnnotation when the expression
    # reads a module score: this process runs AFTER CellTypeAnnotation.
    sobj <- MutateSeuratMeta(sobj, post_mutaters, log = log)
}

# save seurat object
log$info("Saving Seurat object ...")
save_obj(sobj, outfile)
