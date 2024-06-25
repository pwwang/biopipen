source("{{biopipen_dir}}/utils/misc.R")

library(parallel)
library(Seurat)
library(SeuratDisk)
library(rlang)
library(dplyr)

set.seed(8525)

sobjfile = {{in.sobjfile | r}}
outfile = {{out.outfile | r}}
use = {{envs.use | r}}
ident = {{envs.ident | r}}
ref = {{envs.ref | r}}
refnorm = {{envs.refnorm | r}}
ncores = {{envs.ncores | r}}
split_by = {{envs.split_by | r}}
mutaters = {{envs.mutaters | r}}
sctransform_args = {{envs.SCTransform | r: todot="-"}}
normalizedata_args = {{envs.NormalizeData | r: todot="-"}}
findtransferanchors_args = {{envs.FindTransferAnchors | r: todot="-"}}
mappingscore_args = {{envs.MappingScore | r: todot="-"}}
mapquery_args = {{envs.MapQuery | r: todot="-"}}

# See if we have a reference
if (is.null(ref)) {
    stop("No reference provided (envs.ref)")
}

if (is.null(use)) {
    stop("No use provided (envs.use), don't know which column to transfer as cluster")
}

if (is.null(mapquery_args$refdata) || length(mapquery_args$refdata) == 0) {
    mapquery_args$refdata = list()
}

mapquery_args$refdata[[use]] = use

outdir = dirname(outfile)
if (is.null(split_by)) {
    options(future.globals.maxSize = 80000 * 1024^2)
    future::plan(strategy = "multicore", workers = ncores)
}

.expand_dims = function(args, name = "dims") {
    # Expand dims from 30 to 1:30
    if (is.numeric(args[[name]]) && length(args[[name]] == 1)) {
        args[[name]] = 1:args[[name]]
    }
    args
}
findtransferanchors_args = .expand_dims(findtransferanchors_args)

# Load reference
log_info("- Loading reference")
if (endsWith(ref, ".rds") || endsWith(ref, ".RDS")) {
    reference = readRDS(ref)
} else if (endsWith(ref, ".h5ad") || endsWith(ref, ".H5AD")) {
    reference = ReadH5AD(ref)
} else {
    reference = LoadH5Seurat(ref)
}

# check if refdata exists in the reference
for (rname in names(mapquery_args$refdata)) {
    use_name <- mapquery_args$refdata[[rname]]
    # transferring an assay
    if (use_name %in% names(reference)) { next }
    # transferring a metadata column
    if (!use_name %in% colnames(reference@meta.data)) {
        stop(paste0(
            "The reference does not have the column '",
            use_name,
            "' in either assays or metadata. "
        ))
        if (startsWith(use_name, "predicted.")) {
            stop(paste0(
                "Do you mean: ", substring(use_name, 11),
            ))
        }
    }
}

if (refnorm == "auto" && DefaultAssay(reference) == "SCT") {
    refnorm = "SCTransform"
}
log_info("  Normalization method used: {refnorm}")
if (refnorm == "SCTransform") {
    findtransferanchors_args$normalization.method = "SCT"
} else if (refnorm == "NormalizeData") {
    findtransferanchors_args$normalization.method = "LogNormalize"
} else {
    stop("Unknown normalization method: {refnorm}")
}

# Load Seurat object
log_info("- Loading Seurat object")
sobj = readRDS(sobjfile)

if (!is.null(mutaters) && length(mutaters) > 0) {
    log_info("- Applying mutaters")
    sobj@meta.data <- sobj@meta.data %>% mutate(!!!lapply(mutaters, parse_expr))
}

if (!is.null(split_by)) {
    # check if each split has more than 100 cells
    cellno = table(sobj@meta.data[[split_by]])
    cellno = cellno[cellno < 100]
    if (length(cellno) > 0) {
        # stop and print the splits with # cells
        stop(paste0(
            "The following splits have less than 100 cells: \n",
            paste0("- ", names(cellno), ": ", cellno, collapse = "\n"),
            "\n\n",
            "You can use `envs.mutaters` to merge these splits and use `newsplit` as `envs.split_by`: \n",
            "> mutaters = {\n",
            ">   newsplit = \"if_else(oldsplit %in% c('split1', 'split2'), 'mergedsplit', oldsplit)\"\n",
            "> }\n"
        ))
    }
    sobj = SplitObject(sobj, split.by = split_by)
}

# Normalize data
log_info("- Normalizing data")
if (refnorm == "SCTransform") {
    log_info("  Using SCTransform normalization")
    sctransform_args$residual.features = rownames(x = reference)
    if (is.null(split_by)) {
        sctransform_args$object = sobj
        query = do_call(SCTransform, sctransform_args)
    } else {
        query = mclapply(
            X = sobj,
            FUN = function(x) {
                sctransform_args$object = x
                do_call(SCTransform, sctransform_args)
            },
            mc.cores = ncores
        )
        if (any(unlist(lapply(query, class)) == "try-error")) {
            stop(paste0("\nmclapply (SCTransform) error:", query))
        }
    }
} else {
    log_info("  Using NormalizeData normalization")
    if (is.null(split_by)) {
        normalizedata_args$object = sobj
        query = do_call(NormalizeData, normalizedata_args)
    } else {
        query = mclapply(
            X = sobj,
            FUN = function(x) {
                normalizedata_args$object = x
                do_call(NormalizeData, normalizedata_args)
            },
            mc.cores = ncores
        )
        if (any(unlist(lapply(query, class)) == "try-error")) {
            stop(paste0("\nmclapply (NormalizeData) error:", query))
        }
    }
}

# Find anchors between query and reference
log_info("- Finding anchors")
findtransferanchors_args$reference = reference
if (is.null(split_by)) {
    findtransferanchors_args$query = query
    anchors = do_call(FindTransferAnchors, findtransferanchors_args)
} else {
    anchors = mclapply(
        X = query,
        FUN = function(x) {
            findtransferanchors_args$query = x
            do_call(FindTransferAnchors, findtransferanchors_args)
        },
        mc.cores = ncores
    )
    if (any(unlist(lapply(anchors, class)) == "try-error")) {
        stop(paste0("\nmclapply (FindTransferAnchors) error:", anchors))
    }
}

# Map query to reference
log_info("- Mapping query to reference")
mapquery_args$reference = reference
if (is.null(split_by)) {
    mapquery_args$query = query
    mapquery_args$anchorset = anchors
    query = do_call(MapQuery, mapquery_args)
} else {
    query = mclapply(
        X = seq_along(query),
        FUN = function(i) {
            mapquery_args$query = query[[i]]
            mapquery_args$anchorset = anchors[[i]]
            do_call(MapQuery, mapquery_args)
        },
        mc.cores = ncores
    )
    if (any(unlist(lapply(query, class)) == "try-error")) {
        stop(paste0("\nmclapply (MapQuery) error:", query))
    }
}

# Calculating mapping score
log_info("- Calculating mapping score")
mappingscore_sob_msg = paste0(
    "While calculating mapping score, the following error was encountered: \n",
    "subscript out of bounds.  \n\n",
    "You may want to try a smaller `ndim` (default: 50) in `envs.MappingScore`."
)
if (is.null(split_by)) {
    mappingscore_args$anchors = anchors
    mappingscore = tryCatch({
        do_call(MappingScore, mappingscore_args)
    }, error = function(e) {
        if (e$message == "subscript out of bounds") stop(mappingscore_sob_msg)
        stop(e)
    })
} else {
    mappingscore = mclapply(
        X = seq_along(query),
        FUN = function(i) {
            mappingscore_args$anchors = anchors[[i]]
            tryCatch({
                do_call(MappingScore, mappingscore_args)
            }, error = function(e) {
                if (e$message == "subscript out of bounds") stop(mappingscore_sob_msg)
                stop(e)
            })
        },
        mc.cores = ncores
    )
    if (any(unlist(lapply(mappingscore, class)) == "try-error")) {
        stop(paste0("\nmclapply (MappingScore) error:", mappingscore))
    }
}

# Calculate mapping score and add to metadata
log_info("- Adding mapping score to metadata")
if (is.null(split_by)) {
    query = AddMetaData(
        object = query,
        metadata = mappingscore,
        col.name = "mapping.score"
    )
} else {
    query = mclapply(
        X = seq_along(query),
        FUN = function(i) {
            AddMetaData(
                object = query[[i]],
                metadata = mappingscore[[i]],
                col.name = "mapping.score"
            )
        },
        mc.cores = ncores
    )
    if (any(unlist(lapply(query, class)) == "try-error")) {
        stop(paste0("\nmclapply (AddMetaData) error:", query))
    }

    # Combine the results
    log_info("- Merging the results")
    query = merge(query[[1]], query[2:length(query)], merge.dr = "ref.umap")
}

# Add the alias to the metadata for the clusters
log_info("- Adding ident to metadata and set as ident")
query@meta.data = query@meta.data %>% mutate(
    !!sym(ident) := as.factor(!!parse_expr(paste0("predicted.", use)))
)
Idents(query) = ident

# Save
log_info("- Saving result ...")
saveRDS(query, file = outfile)


# ############################
# Some plots
# ############################

# # Plot the UMAP
log_info("- Plotting for transferred data ...")
ref.reduction = mapquery_args$reduction.model %||% "wnn.umap"
for (qname in names(mapquery_args$refdata)) {
    rname <- mapquery_args$refdata[[qname]]

    if (grepl("Array", class(reference[[rname]])) && grepl("Array", class(query[[qname]]))) {
        log_warn("  Skipping transferred array: {qname} -> {rname}")
        next
    }

    log_info("  Plotting transferred data: {qname} -> {rname}")

    ref_p <- DimPlot(
        object = reference,
        reduction = ref.reduction,
        group.by = rname,
        label = TRUE,
        label.size = 3,
        repel = TRUE,
    ) + NoLegend()

    query_p <- DimPlot(
        object = query,
        reduction = "ref.umap",
        group.by = paste0("predicted.", qname),
        label = TRUE,
        label.size = 3,
        repel = TRUE,
    ) + NoLegend()

    png(file.path(outdir, paste0("UMAPs.png")), width = 1400, height = 700, res = 100)
    print(ref_p | query_p)
    dev.off()
}
