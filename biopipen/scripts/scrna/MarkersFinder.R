source("{{biopipen_dir}}/utils/misc.R")
source("{{biopipen_dir}}/utils/mutate_helpers.R")

library(rlang)
library(dplyr)
library(tidyr)
library(tibble)
library(Seurat)
library(enrichR)
library(ggplot2)
library(ggprism)
library(ggrepel)
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
volcano_genes <- {{ envs.volcano_genes | r }}
rest <- {{ envs.rest | r: todot="-" }}
cases <- {{ envs.cases | r: todot="-" }}

if (is.character(volcano_genes) && length(volcano_genes) == 1) {
    volcano_genes <- trimws(strsplit(volcano_genes, ",")[[1]])
}

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
            volcano_genes = volcano_genes,
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
            volcano_genes = volcano_genes,
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

plot_volcano = function(markers, volfile, sig, volgenes) {
    # markers
    #                  gene        p_val avg_log2FC pct.1 pct.2    p_val_adj
    # 1            CCL5 1.883596e-11 -4.8282535 0.359 0.927 4.332270e-09
    # 2        HLA-DQB1 3.667713e-09  6.1543174 0.718 0.098 8.435740e-07
    # 3        HLA-DRB5 1.242993e-07  3.9032231 0.744 0.195 2.858885e-05
    # 4           CD79B 2.036731e-07  4.2748835 0.692 0.146 4.684482e-05
    markers = markers %>%
        mutate(
            Significant = if_else(
                !!parse_expr(sig),
                if_else(avg_log2FC > 0, "Up", "Down"),
                "No"
            ),
            Label = if_else(
                Significant != "No" & (isTRUE(volgenes) | (gene %in% volgenes)),
                gene,
                ""
            )
        )

    p_vol = ggplot(markers, aes(x = avg_log2FC, y = -log10(p_val_adj))) +
        geom_point(aes(color = Significant), alpha = 0.75) +
        scale_color_manual(
            values = c(Up = "#FF3333", Down = "#3333FF", No = "#AAAAAA"),
            labels = c(Up = "Up", Down = "Down", No = "Non-Significant")
        ) +
        geom_text_repel(
            aes(label = Label),
            size = 3,
            color = "#000000",
            box.padding = unit(0.35, "lines"),
            point.padding = unit(0.5, "lines"),
            segment.color = "#000000"
        ) +
        theme_prism() +
        theme(legend.title=element_blank()) +
        labs(
            x = "log2 Fold Change",
            y = "-log10 Adjusted P-value"
        )

    png(volfile, res = 100, height = 800, width = 900)
    print(p_vol)
    dev.off()
}

# Do enrichment analysis for a case using Enrichr
# Args:
#   case: case name
#   markers: markers dataframe
#   sig: The expression to filter significant markers
do_enrich <- function(case, markers, sig, volgenes) {
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
    plot_volcano(markers, file.path(casedir, "volcano.png"), sig, volgenes)
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
    markers <- tryCatch({
        do_call(FindMarkers, args) %>% rownames_to_column("gene")
    }, error = function(e) {
        warning(e$message, immediate. = TRUE)
        data.frame()
    })
    do_enrich(casename, markers, case$sigmarkers, case$volcano_genes)
}

sapply(sort(names(cases)), do_case)
