library(Seurat)
library(rlang)
library(dplyr)
library(tidyselect)
library(biopipen.utils)

srtfile <- {{in.srtobj | r}}
outdir <- {{out.outdir | r}}
joboutdir <- {{job.outdir | r}}
mutaters <- {{ envs.mutaters | r }}
ident <- {{ envs.ident | r }}
group_by <- {{ envs.group_by | default: envs["group-by"] | default: None | r }}  # nolint
each <- {{ envs.each | r }}
dbs <- {{ envs.dbs | r }}
n <- {{ envs.n | r }}
enrich_style <- {{ envs.enrich_style | r }}
sset <- {{ envs.subset | r }}
enrich_plots_defaults <- {{ envs.enrich_plots_defaults | r }}
enrich_plots <- {{ envs.enrich_plots | r }}
cases <- {{ envs.cases | r: todot = "-" }}  # nolint

set.seed(8525)
log <- get_logger()
reporter <- get_reporter()

log$info("Reading Seurat object ...")
srtobj <- read_obj(srtfile)
assay <- DefaultAssay(srtobj)

if (!is.null(mutaters) && length(mutaters) > 0) {
    log$info("Mutating meta data ...")
    srtobj@meta.data <- srtobj@meta.data %>%
        mutate(!!!lapply(mutaters, parse_expr))
}

enrich_plots <- lapply(enrich_plots, function(x) {
    list_update(enrich_plots_defaults, x)
})
defaults <- list(
    ident = ident,
    group_by = group_by,
    each = each,
    dbs = dbs,
    n = n,
    enrich_style = enrich_style,
    enrich_plots = enrich_plots,
    enrich_plots_defaults = enrich_plots_defaults,
    subset = sset
)

cases <- expand_cases(cases, defaults, default_case = "Top Expressing Genes", post = function(name, case) {
    outcases <- list()
    if (is.null(case$each) || is.na(case$each) || nchar(case$each) == 0 || isFALSE(each)) {
        case$enrich_plots <- lapply(
            case$enrich_plots,
            function(x) { list_update(case$enrich_plots_defaults, x) }
        )
        case$enrich_plots_defaults <- NULL

        outcases[[name]] <- case
    } else {
        eachs <- if (!is.null(case$subset)) {
            srtobj@meta.data %>%
                filter(!!parse_expr(case$subset)) %>%
                pull(case$each) %>% na.omit() %>% unique() %>% as.vector()
        } else {
            srtobj@meta.data %>%
                pull(case$each) %>% na.omit() %>% unique() %>% as.vector()
        }

        if (length(cases) == 0 && name == "Top Expressing Genes") {
            name <- case$each
        }

        for (each in eachs) {
            newname <- paste0(name, " - ", each)
            newcase <- case
            newcase$each_name <- case$each
            newcase$each <- each

            if (!is.null(case$subset)) {
                newcase$subset <- paste0(case$subset, " & ", bQuote(case$each), " == '", each, "'")
            } else {
                newcase$subset <- paste0(bQuote(case$each), " == '", each, "'")
            }

            newcase$enrich_plots <- lapply(
                case$enrich_plots,
                function(x) { list_update(case$enrich_plots_defaults, x) }
            )
            newcase$enrich_plots_defaults <- NULL

            outcases[[newname]] <- newcase
        }
    }

    outcases
})

log$info("Running cases ...")

process_markers <- function(markers, info, case) {
    # Save markers
    write.table(markers, file.path(info$prefix, "top_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    reporter$add2(
        list(
            name = "Table",
            contents = list(
                list(kind = "descr", content = "Showing top expressing genes ordered by their expression descendingly."),
                list(kind = "table", src = file.path(info$prefix, "top_genes.tsv"), data = list(nrows = 100))
            )
        ),
        hs = c(info$section, info$name),
        hs2 = paste0("Top Genes"),
        ui = "tabs"
    )

    enrich <- RunEnrichment(
        markers$gene,
        dbs = case$dbs, style = case$enrich_style)

    write.table(enrich, file.path(info$prefix, "enrich.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
    reporter$add2(
        list(
            name = "Table",
            contents = list(list(kind = "table", src = file.path(info$prefix, "enrich.tsv"), data = list(nrows = 100)))
        ),
        hs = c(info$section, info$name),
        hs2 = "Enrichment Analysis",
        ui = "tabs"
    )

    # Visualize enriched terms
    if (length(case$enrich_plots) > 0) {
        for (db in case$dbs) {
            plots <- list()
            for (plotname in names(case$enrich_plots)) {
                plotargs <- case$enrich_plots[[plotname]]
                plotargs$data <- enrich[enrich$Database == db, , drop = FALSE]

                p <- do_call(VizEnrichment, plotargs)

                outprefix <- file.path(info$prefix, paste0("enrich.", slugify(db), ".", slugify(plotname)))
                if (plotargs$plot_type == "bar") {
                    attr(p, "height") <- attr(p, "height") / 1.5
                }
                save_plot(p, outprefix, plotargs$devpars, formats = "png")
                plots[[length(plots) + 1]] <- reporter$image(outprefix, c(), FALSE)
            }
            reporter$add2(
                list(name = db, contents = plots),
                hs = c(info$section, info$name),
                hs2 = "Enrichment Analysis",
                ui = "tabs"
            )
        }
    }
}


run_case <- function(name) {
    log$info("Case: {name} ...")
    case <- cases[[name]]

    log$info("- Subsetting cells and calculating average expression ...")
    if (!is.null(case$subset)) {
        subobj <- filter(srtobj, !!parse_expr(case$subset))
    } else {
        subobj <- srtobj
    }
    case$group_by <- case$group_by %||% GetIdentityColumn(srtobj)
    if (is.null(case$ident)) {
        case$ident <- as.character(unique(subobj@meta.data[[case$group_by]]))
    }
    avgexpr <- AverageExpression(
        subobj,
        group_by = case$group_by,
        assays = assay
    )[[assay]]
    # https://github.com/satijalab/seurat/issues/7893
    colnames(avgexpr) <- as.character(unique(subobj@meta.data[[case$group_by]]))
    avgexpr <- avgexpr[, case$ident, drop = FALSE]

    for (idt in case$ident) {
        log$info("- Processing {idt} ...")
        info <- case_info(paste0(name, "::", idt), outdir, create = TRUE)
        expr <- avgexpr[, idt, drop = FALSE]
        expr <- expr[order(expr, decreasing = TRUE), , drop = FALSE]
        expr <- expr[1:min(case$n, nrow(expr)), , drop = FALSE]
        expr <- as.data.frame(expr)
        expr$gene <- rownames(expr)
        colnames(expr) <- c("avg_expr", "gene")
        expr <- expr[, c("gene", "avg_expr"), drop = FALSE]

        log$info("  Performing enrichment analysis ...")
        process_markers(expr, info, case = list(
            ident = idt,
            dbs = case$dbs,
            enrich_style = case$enrich_style,
            enrich_plots = case$enrich_plots
        ))
    }

    invisible()
}

sapply(names(cases), run_case)

reporter$save(joboutdir)
