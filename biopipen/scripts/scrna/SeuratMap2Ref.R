source("{{biopipen_dir}}/utils/misc.R")

library(Seurat)
library(SeuratDisk)
library(rlang)
library(dplyr)

set.seed(8525)

sobjfile = {{in.sobjfile | r}}
outfile = {{out.outfile | r}}
use = {{envs.use | r}}
name = {{envs.name | r}}
ref = {{envs.ref | r}}
ncores = {{envs.ncores | r}}
sctransform_args = {{envs.SCTransform | r: todot="-"}}
findtransferanchors_args = {{envs.FindTransferAnchors | r: todot="-"}}
mappingscore_args = {{envs.MappingScore | r: todot="-"}}
mapquery_args = {{envs.MapQuery | r: todot="-"}}

if (is.null(mapquery_args$refdata) || length(mapquery_args$refdata) == 0) {
    stop("No refdata provided for MapQuery (envs.MapQuery.refdata)")
}

# See if we have a reference
if (is.null(ref)) {
    stop("No reference provided (envs.ref)")
}

if (is.null(use)) {
    if (length(mapquery_args$refdata) == 1) {
        use = paste0("predicted.", names(mapquery_args$refdata)[1])
    } else {
        stop("No use provided (envs.use), don't know which column to use as cluster")
    }
}

outdir = dirname(outfile)
options(future.globals.maxSize = 80000 * 1024^2)
plan(strategy = "multicore", workers = ncores)

.expand_dims = function(args, name = "dims") {
    # Expand dims from 30 to 1:30
    if (is.numeric(args[[name]]) && length(args[[name]] == 1)) {
        args[[name]] = 1:args[[name]]
    }
    args
}
findtransferanchors_args = .expand_dims(findtransferanchors_args)

# Load reference
log_info("Loading reference")
if (endsWith(ref, ".rds") || endsWith(ref, ".RDS")) {
    reference = readRDS(ref)
} else if (endsWith(ref, ".h5ad") || endsWith(ref, ".H5AD")) {
    reference = ReadH5AD(ref)
} else {
    reference = LoadH5Seurat(ref)
}

# Load Seurat object
log_info("Loading Seurat object")
sobj = readRDS(sobjfile)

# Normalize data
log_info("Normalizing data")
sctransform_args$object = sobj
sctransform_args$residual.features = rownames(x = reference)
query = do_call(SCTransform, sctransform_args)

# Find anchors between query and reference
log_info("Finding anchors")
findtransferanchors_args$reference = reference
findtransferanchors_args$query = query
anchors = do_call(FindTransferAnchors, findtransferanchors_args)

# Map query to reference
log_info("Mapping query to reference")
mapquery_args$reference = reference
mapquery_args$query = query
mapquery_args$anchorset = anchors
query = do_call(MapQuery, mapquery_args)

# Calculating mapping score
log_info("Calculating mapping score")
mappingscore_args$anchors = anchors
mappingscore = tryCatch({
    do_call(MappingScore, mappingscore_args)
}, error = function(e) {
    if (e$message == "subscript out of bounds") {
        stop(paste0(
            "While calculating mapping score, the following error was encountered: \n",
            "subscript out of bounds.  \n\n",
            "You may want to try a smaller `ndim` (default: 50) in `envs.MappingScore`."
        ))
    }
    stop(e)
})

# Calculate mapping score and add to metadata
log_info("Calculating mapping score")
query = AddMetaData(
  object = query,
  metadata = mappingscore,
  col.name = "mapping.score"
)

# Add the alias to the metadata for the clusters
log_info("Adding name to metadata and set as ident")
query@meta.data = query@meta.data %>% mutate(!!sym(name) := as.factor(!!parse_expr(use)))
Idents(query) = name

# Save
log_info("Saving")
saveRDS(query, file = outfile)


# ############################
# Some plots
# ############################

# # Plot the UMAP
log_info("Plotting")
for (refname in names(mapquery_args$refdata)) {
    if (refname == "predicted_ADT") {
        next
    }
    reduction = if (is.null(mapquery_args$reduction.model)) "wnn.umap" else mapquery_args$reduction.model
    print(paste0("Plotting UMAP for ", refname))
    p = DimPlot(
        object = reference,
        reduction = reduction,
        group.by = refname,
        label = TRUE,
        label.size = 3,
        repel = TRUE,
    ) + NoLegend()

    png(file.path(outdir, paste0("Reference_UMAP_", refname, ".png")), width = 1000, height = 1000, res = 100)
    print(p)
    dev.off()

    p = DimPlot(
        object = query,
        reduction = "ref.umap",
        group.by = paste0("predicted.", refname),
        label = TRUE,
        label.size = 3,
        repel = TRUE,
    ) + NoLegend()

    png(file.path(outdir, paste0("Query_UMAP_", refname, ".png")), width = 1000, height = 1000, res = 100)
    print(p)
    dev.off()
}
