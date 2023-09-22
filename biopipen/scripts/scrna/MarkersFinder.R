source("{{biopipen_dir}}/utils/misc.R")
source("{{biopipen_dir}}/utils/mutate_helpers.R")

library(rlang)
library(dplyr)
library(tidyr)
library(tibble)
library(Seurat)
library(enrichR)
library(ggplot2)
library(future)
library(tidyseurat)

setEnrichrSite("Enrichr")

srtfile <- {{ in.srtobj | quote }}
outdir <- {{ out.outdir | quote }}
ncores <- {{ envs.ncores | int }}
mutaters <- {{ envs.mutaters | r }}
ident.1 <- {{ envs["ident-1"] | r }}
ident.2 <- {{ envs["ident-2"] | r }}
group.by <- {{ envs["group-by"] | r }}
each <- {{ envs.each | r }}
prefix_each <- {{ envs.prefix_each | r }}
section <- {{ envs.section | r }}
dbs <- {{ envs.dbs | r }}
sigmarkers <- {{ envs.sigmarkers | r }}
rest <- {{ envs.rest | r: todot="-" }}
cases <- {{ envs.cases | r: todot="-" }}

set.seed(8525)
if (ncores > 1) {
    options(future.globals.maxSize = 80000 * 1024^2)
    plan(strategy = "multicore", workers = ncores)
}

print("- Reading Seurat object ...")
srtobj <- readRDS(srtfile)

print("- Mutate meta data if needed ...")
if (!is.null(mutaters) && length(mutaters) > 0) {
    srtobj@meta.data <- srtobj@meta.data %>%
        mutate(!!!lapply(mutaters, parse_expr))
}

print("- Expanding cases ...")
if (is.null(cases) || length(cases) == 0) {
    cases <- list(
        DEFAULT = list(
            ident.1 = ident.1,
            ident.2 = ident.2,
            group.by = group.by,
            each = each,
            prefix_each = prefix_each,
            section = section,
            dbs = dbs,
            sigmarkers = sigmarkers,
            rest = rest
        )
    )
} else {
    for (name in names(cases)) {
        case <- list_setdefault(
            cases[[name]],
            ident.1 = ident.1,
            ident.2 = ident.2,
            group.by = group.by,
            each = each,
            prefix_each = prefix_each,
            section = section,
            dbs = dbs,
            sigmarkers = sigmarkers,
            rest = rest
        )
        case$rest <- list_setdefault(case$rest, rest)
        cases[[name]] <- case
    }
}
# Expand each and with ident.1
#  list(Cluster0 = list(each = "Sample", group.by = "seurat_clusters", ident.1 = "0"))
# to
#  list(
#   `Sample-Sample1:Cluster0` = list(...),
#   `Sample-Sample2:Cluster0` = list(...),
#   ...
#  )
# Expand each and without ident.1
#  list(Cluster = list(each = "Sample", group.by = "seurat_clusters"))
# to
#  list(
#   `Sample-Sample1-Cluster:0` = list(...),
#   `Sample-Sample1-Cluster:1` = list(...),
#   ...
#   `Sample2-Cluster:0` = list(...),
#   `Sample2-Cluster:1` = list(...),
#   ...
#  )
# If no each, and not ident.1
#  list(Cluster = list(group.by = "seurat_clusters"))
# to
#  list(
#   `Cluster:0` = list(...),
#   `Cluster:1` = list(...),
#   ...
#  )
# Otherwise if section is specified, the case name will be changed to `section:case`

newcases <- list()
for (name in names(cases)) {
    case <- cases[[name]]
    if (is.null(case$each) && !is.null(case$ident.1)) {
        newcases[[paste0(case$section, ":", name)]] <- case
    } else if (is.null(case$each)) {
        # is.null(case$ident.1)
        idents <- srtobj@meta.data %>% pull(case$group.by) %>% unique() %>% na.omit()
        for (ident in idents) {
            newcases[[paste0(name, ":", ident)]] <- case
            newcases[[paste0(name, ":", ident)]]$ident.1 <- ident
        }
    } else {
        eachs <- srtobj@meta.data %>% pull(case$each) %>% unique() %>% na.omit()
        for (each in eachs) {
            by = make.names(paste0(".", name, "_", each))
            srtobj@meta.data = srtobj@meta.data %>% mutate(
                !!sym(by) := if_else(
                    !!sym(case$each) == each,
                    !!sym(case$group.by),
                    NA
                )
            )
            if (is.null(case$ident.1)) {
                idents <- srtobj@meta.data %>% pull(case$group.by) %>% unique() %>% na.omit()
                for (ident in idents) {
                    kname <- if (name == "DEFAULT") "" else paste0("-", name)
                    key <- paste0(each, kname, ":", ident)
                    if (case$prefix_each) {
                        key <- paste0(case$each, "-", key)
                    }
                    newcases[[key]] <- case
                    newcases[[key]]$ident.1 <- ident
                    newcases[[key]]$group.by <- by
                }
            } else {
                key <- paste0(case$each, ":", each)
                if (name != "DEFAULT") {
                    key <- paste0(key, " - ", name)
                }
                newcases[[key]] <- case
                newcases[[key]]$group.by <- by
            }
        }
    }
}
cases <- newcases


# Do enrichment analysis for a case using Enrichr
# Args:
#   case: case name
#   markers: markers dataframe
#   sig: The expression to filter significant markers
do_enrich <- function(case, markers, sig) {
    print(paste("  Running enrichment for case:", case))
    parts <- strsplit(case, ":")[[1]]
    sec <- parts[1]
    case <- paste0(parts[-1], collapse = ":")
    casedir <- file.path(outdir, sec, case)
    dir.create(casedir, showWarnings = FALSE, recursive = TRUE)
    if (nrow(markers) == 0) {
        print(paste("  No markers found for case:", case))
        cat("No markers found.", file = file.path(casedir, "error.txt"))
        return()
    }
    markers_sig <- markers %>% filter(!!parse_expr(sig))
    if (nrow(markers_sig) == 0) {
        print(paste("  No significant markers found for case:", case))
        cat("No significant markers.", file = file.path(casedir, "error.txt"))
        return()
    }
    write.table(
        markers_sig,
        file.path(casedir, "markers.txt"),
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE
    )
    if (nrow(markers_sig) < 5) {
        for (db in dbs) {
            write.table(
                data.frame(Warning = "Not enough significant markers."),
                file.path(casedir, paste0("Enrichr-", db, ".txt")),
                sep = "\t",
                row.names = FALSE,
                col.names = TRUE,
                quote = FALSE
            )
            png(
                file.path(casedir, paste0("Enrichr-", db, ".png")),
                res = 100, height = 200, width = 1000
            )
            print(
                ggplot() +
                    annotate(
                        "text",
                        x = 1,
                        y = 1,
                        label = "Not enough significant markers."
                    ) +
                    theme_classic()
            )
            dev.off()
        }
    } else {
        enriched <- enrichr(markers_sig$gene, dbs)
        for (db in dbs) {
            write.table(
                enriched[[db]],
                file.path(casedir, paste0("Enrichr-", db, ".txt")),
                sep = "\t",
                row.names = FALSE,
                col.names = TRUE,
                quote = FALSE
            )
            png(
                file.path(casedir, paste0("Enrichr-", db, ".png")),
                res = 100, height = 1000, width = 1000
            )
            print(plotEnrich(enriched[[db]], showTerms = 20, title = db))
            dev.off()
        }
    }
}


do_case <- function(casename) {
    cat(paste("- Dealing with case:", casename, "...\n"))
    case <- cases[[casename]]
    # ident1
    # ident2
    # groupby
    # each  # expanded
    # prefix_each
    # dbs
    # sigmarkers
    # rest
    args <- case$rest
    args$group.by <- case$group.by
    args$ident.1 <- case$ident.1
    args$ident.2 <- case$ident.2
    idents <- srtobj@meta.data %>% pull(case$group.by) %>% unique()
    if (anyNA(idents)) {
        args$object <- srtobj %>% filter(!is.na(!!sym(case$group.by)))
    } else {
        args$object <- srtobj
    }
    markers <- do_call(FindMarkers, args) %>% rownames_to_column("gene")
    do_enrich(casename, markers, case$sigmarkers)
}

sapply(sort(names(cases)), do_case)
