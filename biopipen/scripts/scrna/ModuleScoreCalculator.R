source("{{biopipen_dir}}/utils/misc.R")
library(Seurat)
library(dplyr)

sobjfile <- {{in.srtobj | r}}
outfile <- {{out.rdsfile | r}}
defaults <- {{envs.defaults | r}}
modules <- {{envs.modules | r}}

# load seurat object
log_info("Loading Seurat object ...")
sobj <- readRDS(sobjfile)

aggs <- list(
    mean = mean,
    median = median,
    sum = sum,
    max = max,
    min = min,
    sd = sd,
    var = var
)

for (key in names(modules)) {
    if (is.null(modules[[key]])) {
        modules[[key]] <- list()
    }

    module <- list_update(defaults, modules[[key]])
    if (is.null(module$features) || length(module$features) == 0) {
        stop(paste0("Module '", key, "' has no features"))
    }

    keep <- module$keep
    agg <- aggs[[module$agg]]
    module$keep <- NULL
    module$agg <- NULL
    log_info("Calculating module '{key}' ...")
    is_cc <- FALSE
    if (!is.null(module$kind) && module$kind %in% c("diffmap", "diffusion_map")) {
        library(destiny)
        features <- module$features
        if (is.null(features)) { features <- 2 }
        if (is.null(module$verbose)) { module$verbose <- TRUE }
        module$features <- NULL
        module$kind <- NULL

        if (!is.null(module$n_pcs)) {
            log_info("- Using cell embeddings from PCA reduction ...")
            module$data <- Embeddings(sobj, reduction = "pca")
            if (module$n_pcs > ncol(module$data)) {
                log_warn("- `n_pcs` ({module$n_pcs}) is larger than the number of PCs, using all {ncol(module$data)} PCs ...")
            }
            module$data <- module$data[, 1:min(module$n_pcs, ncol(module$data))]
            module$n_pcs <- NULL
        } else {
            log_info("- Using assay data ...")
            module$data <- GetAssayData(sobj, layer = "data")
        }

        log_info("- Calculating diffusion map ...")
        dm <- do_call(DiffusionMap, module)
        ev <- eigenvectors(dm)

        log_info("- Creating DimReduc object ...")
        sobj[[key]] <- CreateDimReducObject(
            embeddings = data.matrix(as.data.frame(ev[, 1:features])),
            key = paste0(key, "_")
        )

        # add to meta.data
        log_info("- Adding to meta.data ...")
        sobj <- AddMetaData(
            sobj,
            sobj[[key]]@cell.embeddings,
            col.name = colnames(sobj[[key]]@cell.embeddings)
        )

        next
    }

    module$object <- sobj
    if (length(module$features) == 1 && module$features == "cc.genes") {
        is_cc <- TRUE
        module$features <- NULL
        module$s.features <- cc.genes$s.genes
        module$g2m.features <- cc.genes$g2m.genes
    } else if (length(module$features) == 1 && module$features == "cc.genes.updated.2019") {
        is_cc <- TRUE
        module$features <- NULL
        module$s.features <- cc.genes.updated.2019$s.genes
        module$g2m.features <- cc.genes.updated.2019$g2m.genes
    } else {
        module$name <- key
        if (length(module$features) == 1) {
            module$features <- trimws(strsplit(module$features, ",")[[1]])
        }
        module$features <- list(module$features)
    }
    if (isTRUE(is_cc)) {
        sobj <- do_call(CellCycleScoring, module)
    } else {
        tryCatch({
            sobj <- do_call(AddModuleScore, module)
        }, error = function(e) {
            if (grepl("cannot take a sample larger than", e$message)) {
                stop(paste0(
                    "Module '", key, "': ",
                    e$message,
                    " (try a smaller `ctrl`?)"
                ))
            } else {
                stop(e)
            }
        })

        sobj[[key]] <- sobj@meta.data %>%
            rowwise() %>%
            mutate(
                !!sym(key) := agg(
                    c_across(matches(paste0("^", key, "\\d+$"))),
                    na.rm = TRUE
                )
            ) %>%
            ungroup() %>%
            pull(key)

        if (!isTRUE(keep)) {
            sobj@meta.data <- sobj@meta.data %>%
                select(-matches(paste0("^", key, "\\d+$")))
        }
    }
}

# save seurat object
print("Saving Seurat object ...")
saveRDS(sobj, outfile)
