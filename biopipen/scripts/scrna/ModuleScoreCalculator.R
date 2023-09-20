source("{{biopipen_dir}}/utils/misc.R")
library(Seurat)
library(dplyr)

sobjfile <- {{in.srtobj | r}}
outfile <- {{out.rdsfile | r}}
defaults <- {{envs.defaults | r}}
modules <- {{envs.modules | r}}

# load seurat object
print("Loading Seurat object ...")
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
    module$object <- sobj
    if (is.null(module$features) || length(module$features) == 0) {
        stop(paste0("Module '", key, "' has no features"))
    }

    keep <- module$keep
    agg <- aggs[[module$agg]]
    module$keep <- NULL
    module$agg <- NULL
    print(paste0("Calculating module '", key, "' ..."))
    is_cc <- FALSE
    if (module$features == "cc.genes") {
        is_cc <- TRUE
        module$features <- NULL
        module$s.features <- cc.genes$s.genes
        module$g2m.features <- cc.genes$g2m.genes
    } else if (module$features == "cc.genes.updated.2019") {
        is_cc <- TRUE
        module$features <- NULL
        module$s.features <- cc.genes.updated.2019$s.genes
        module$g2m.features <- cc.genes.updated.2019$g2m.genes
    } else {
        module$name <- key
        module$features <- trimws(strsplit(module$features, ",")[[1]])
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
