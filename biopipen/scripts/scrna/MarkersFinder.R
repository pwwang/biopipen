source("{{biopipen_dir}}/utils/misc.R")
source("{{biopipen_dir}}/utils/caching.R")
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
library(ggVennDiagram)
library(UpSetR)

log_info("Setting up EnrichR ...")
setEnrichrSite("Enrichr")

srtfile <- {{ in.srtobj | quote }}
outdir <- {{ out.outdir | quote }}
joboutdir <- {{ job.outdir | quote }}
ncores <- {{ envs.ncores | int }}
mutaters <- {{ envs.mutaters | r }}
ident.1 <- {{ envs["ident-1"] | r }}
ident.2 <- {{ envs["ident-2"] | r }}
group.by <- {{ envs["group-by"] | r }}
each <- {{ envs.each | r }}
prefix_each <- {{ envs.prefix_each | r }}
prefix_group <- {{ envs.prefix_group | r }}
section <- {{ envs.section | r }}
dbs <- {{ envs.dbs | r }}
assay <- {{ envs.assay | r }}
sigmarkers <- {{ envs.sigmarkers | r }}
volcano_genes <- {{ envs.volcano_genes | r }}
subset <- {{ envs.subset | r }}
rest <- {{ envs.rest | r: todot="-" }}
dotplot <- {{ envs.dotplot | r: todot="-" }}
cases <- {{ envs.cases | r: todot="-", skip=1 }}
overlapping_defaults <- {{ envs.overlap_defaults | r }}
overlapping <- {{ envs.overlap | r }}
cache <- {{ envs.cache | r }}

if (isTRUE(cache)) { cache <- joboutdir }

# expand overlapping
for (sec in names(overlapping)) {
    overlapping[[sec]] <- list_update(overlapping_defaults, overlapping[[sec]])
}
overlapping_sections <- names(overlapping)

overlaps <- list()
if (is.character(volcano_genes) && length(volcano_genes) == 1) {
    volcano_genes <- trimws(strsplit(volcano_genes, ",")[[1]])
}

set.seed(8525)
if (ncores > 1) {
    options(future.globals.maxSize = 80000 * 1024^2)
    plan(strategy = "multicore", workers = ncores)
}

log_info("- Reading Seurat object ...")
srtobj <- readRDS(srtfile)
defassay <- DefaultAssay(srtobj)
if (defassay == "SCT" && !"PrepSCTFindMarkers" %in% names(srtobj@commands)) {
    log_warn("  SCTransform used but PrepSCTFindMarkers not applied, running ...")

    srtobj <- PrepSCTFindMarkers(srtobj)
    # compose a new SeuratCommand to record it to srtobj@commands
    commands <- names(srtobj@commands)
    scommand <- srtobj@commands[[commands[length(commands)]]]
    scommand@name <- "PrepSCTFindMarkers"
    scommand@time.stamp <- Sys.time()
    scommand@assay.used <- "SCT"
    scommand@call.string <- "PrepSCTFindMarkers(object = srtobj)"
    scommand@params <- list()
    srtobj@commands$PrepSCTFindMarkers <- scommand
}

if (!is.null(mutaters) && length(mutaters) > 0) {
    log_info("- Mutating meta data ...")
    srtobj@meta.data <- srtobj@meta.data %>%
        mutate(!!!lapply(mutaters, parse_expr))
}

defaults <- list(
    ident.1 = ident.1,
    ident.2 = ident.2,
    group.by = group.by,
    each = each,
    prefix_each = prefix_each,
    prefix_group = prefix_group,
    section = section,
    dbs = dbs,
    assay = assay %||% defassay,
    subset = subset,
    sigmarkers = sigmarkers,
    volcano_genes = volcano_genes,
    dotplot = dotplot,
    rest = rest
)

expand_each <- function(name, case) {
    outcases <- list()
    no_each <- is.null(case$each) || nchar(case$each) == 0
    if (no_each && !is.null(case$ident.1)) {
        # single cases, no need to expand
        if (is.null(case$section) || case$section == "DEFAULT") {
            outcases[[name]] <- case
        } else {
            outcases[[paste0(case$section, "::", name)]] <- case
        }
    } else {  # !no_each || is.null(case$ident.1)
        if (!is.null(case$section) && case$section != "DEFAULT") {
            log_warn("  Ignoring `section` in case `{name}` that will be expanded (`each` is set or `ident-1` is not set).")
            case$section <- NULL
        }
        if (no_each) {  # is.null(ident.1)
            # no each and no ident.1, use FindAllMarkers
            key <- paste0(name, "::", name)
            outcases[[key]] <- case
            outcases[[key]]$section <- name
            outcases[[key]]$findall <- TRUE
        } else if (!no_each) {
            # expand each
            if (is.null(case$subset)) {
                eachs <- srtobj@meta.data %>%
                    pull(case$each) %>% na.omit() %>% unique() %>% as.vector()
            } else {
                eachs <- srtobj@meta.data %>% dplyr::filter(!!parse_expr(case$subset)) %>%
                    pull(case$each) %>% na.omit() %>% unique() %>% as.vector()
            }
            for (each in eachs) {
                by <- make.names(paste0("..", name, "_", case$each,"_", each))
                srtobj@meta.data <<- srtobj@meta.data %>% mutate(
                    !!sym(by) := if_else(
                        !!sym(case$each) == each,
                        !!sym(case$group.by),
                        NA
                    )
                )
                if (isTRUE(case$prefix_each)) {
                    key <- paste0(name, "::", case$each, " - ", each)
                } else {
                    key <- paste0(name, "::", each)
                }
                outcases[[key]] <- case
                outcases[[key]]$section <- name
                outcases[[key]]$group.by <- by
                if (is.null(case$ident.1)) {
                    outcases[[key]]$findall <- TRUE
                }
            }
        }
    }
    outcases
}

log_info("- Expanding cases ...")
cases <- expand_cases(cases, defaults, expand_each)

plot_volcano = function(markers, volfile, sig, volgenes) {
    # markers
    #                  gene        p_val avg_log2FC pct.1 pct.2    p_val_adj
    # 1            CCL5 1.883596e-11 -4.8282535 0.359 0.927 4.332270e-09
    # 2        HLA-DQB1 3.667713e-09  6.1543174 0.718 0.098 8.435740e-07
    # 3        HLA-DRB5 1.242993e-07  3.9032231 0.744 0.195 2.858885e-05
    # 4           CD79B 2.036731e-07  4.2748835 0.692 0.146 4.684482e-05
    log_info("  Plotting volcano plot ...")
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
        theme(legend.title=element_blank(), plot.margin=unit(c(1,1,1,1), "cm")) +
        labs(
            x = "log2 Fold Change",
            y = "-log10 Adjusted P-value"
        )

    png(volfile, res = 100, height = 1200, width = 900)
    print(p_vol)
    dev.off()
}

# Do enrichment analysis for a case using Enrichr
# Args:
#   case: case name
#   markers: markers dataframe
#   sig: The expression to filter significant markers
do_enrich <- function(info, markers, sig, volgenes) {
    log_info("  Running enrichment for case: {info$casename}")

    if (nrow(markers) == 0) {
        log_warn("  No markers found for case: {info$casename}")
        return(NULL)
    }

    plot_volcano(markers, file.path(info$casedir, "volcano.png"), sig, volgenes)

    markers_sig <- markers %>% filter(!!parse_expr(sig)) %>% arrange(p_val_adj)
    if (nrow(markers_sig) == 0) {
        log_warn("  No significant markers found.")
        return(NULL)
    }

    write.table(
        markers_sig,
        file.path(info$casedir, "markers.txt"),
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE
    )
    if (nrow(markers_sig) < 5) {
        log_warn("  Too few significant markers found for case: {info$casename}")
    } else {
        enriched <- enrichr(unique(markers_sig$gene), dbs)
        for (db in dbs) {
            write.table(
                enriched[[db]],
                file.path(info$casedir, paste0("Enrichr-", db, ".txt")),
                sep = "\t",
                row.names = FALSE,
                col.names = TRUE,
                quote = FALSE
            )
            if (nrow(enriched[[db]]) == 0) {
                log_warn("  No enrichment found for case: {info$casename} - {db}")
                next
            }
            png(
                file.path(info$casedir, paste0("Enrichr-", db, ".png")),
                res = 100, height = 1000, width = 1000
            )
            print(
                plotEnrich(enriched[[db]], showTerms = 20, title = db) +
                theme_prism()
            )
            dev.off()
        }
    }
    unique(markers_sig$gene)
}

do_dotplot <- function(info, siggenes, dotplot, args) {
    max_dotplot_features <- dotplot$maxgenes %||% 20
    dotplot$maxgenes <- NULL
    if (length(siggenes) > max_dotplot_features) {
        log_debug("  Too many significant markers ({length(siggenes)}), using first {max_dotplot_features} for dotplot")
        siggenes <- siggenes[1:max_dotplot_features]
    }
    dotplot_devpars <- dotplot$devpars
    dotplot$devpars <- NULL
    dotplot$object <- args$object
    dotplot$features <- siggenes
    dotplot$group.by <- args$group.by
    dotplot_width <- dotplot_devpars$width %||%
        ifelse(length(siggenes) <= 20, length(siggenes) * 60, min(1000, length(siggenes)) * 30)
    dotplot_height <- dotplot_devpars$height %||% 600
    dotplot_res <- dotplot_devpars$res %||% 100
    dotplot_file <- file.path(info$casedir, "dotplot.png")
    png(dotplot_file, res = dotplot_res, width = dotplot_height, height = dotplot_width)
    # rotate x axis labels
    print(
        do_call(DotPlot, dotplot) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
        coord_flip()
    )
    dev.off()
}

add_case_report <- function(info, sigmarkers, siggenes) {
    h1 = info$h1
    h2 = info$h2
    if (is.null(siggenes) || length(siggenes) == 0) {
        add_report(
            list(
                kind = "error",
                content = "No significant markers found."
            ),
            h1 = h1,
            h2 = ifelse(h2 == "#", "Markers", h2),
            h3 = ifelse(h2 == "#", "#", "Markers"),
            ui = "flat"
        )
    } else {
        add_report(
            list(
                title = "Significant Markers",
                ui = "flat",
                contents = list(
                    list(
                        kind = "descr",
                        content = paste0(
                            "The markers are found using Seurat's FindMarkers function, ",
                            "and filtered by: ",
                            html_escape(sigmarkers)
                        )
                    ),
                    list(
                        kind = "table",
                        data = list(nrows = 100),
                        src = file.path(info$casedir, "markers.txt")
                    )
                )
            ),
            list(
                title = "Volcano Plot",
                ui = "flat",
                contents = list(
                    list(
                        kind = "img",
                        src = file.path(info$casedir, "volcano.png")
                    )
                )
            ),
            list(
                title = "Dot Plot",
                ui = "flat",
                contents = list(
                    list(
                        kind = "img",
                        src = file.path(info$casedir, "dotplot.png")
                    )
                )
            ),
            h1 = h1,
            h2 = ifelse(h2 == "#", "Markers", h2),
            h3 = ifelse(h2 == "#", "#", "Markers"),
            ui = "tabs"
        )

        add_report(
            list(
                kind = "descr",
                content = paste0(
                    "The enrichment analysis is done using Enrichr. ",
                    "The significant markers are used as input. "
                )
            ),
            list(
                kind = "enrichr",
                dir = info$casedir
            ),
            h1 = h1,
            h2 = ifelse(h2 == "#", "Enrichment Analysis", h2),
            h3 = ifelse(h2 == "#", "#", "Enrichment Analysis"),
            ui = "flat"
        )
    }
}

ensure_sobj <- function(expr, allow_empty) {
    tryCatch({ expr }, error = function(e) {
        if (allow_empty) {
            log_warn("  Ignoring this case: {e$message}")
            return(NULL)
        } else {
            stop(e)
        }
    })
}

do_case_findall <- function(casename) {
    # casename
    ## Cluster::Cluster
    info <- casename_info(casename, cases, outdir, create = FALSE)
    if (info$section %in% overlapping_sections) {
        stop(paste0("  Can't do overlapping analysis for case without `ident-1` set: ", casename))
    }

    case <- cases[[casename]]
    log_info("  Using FindAllMarkers for case: {casename}...")
    args <- case$rest
    args$assay <- case$assay
    args$group.by <- case$group.by
    # args$logfc.threshold <- args$logfc.threshold %||% 0
    # args$min.cells.group <- args$min.cells.group %||% 1
    # args$min.cells.feature <- args$min.cells.feature %||% 1
    # args$min.pct <- args$min.pct %||% 0
    allow_empty = startsWith(case$group.by, "..")
    if (!is.null(case$subset)) {
        args$object <- ensure_sobj({
            srtobj %>% filter(!!parse_expr(case$subset) & !is.na(!!sym(case$group.by)))
        }, allow_empty)
        if (is.null(args$object)) { return() }
    } else {
        args$object <- ensure_sobj({
            srtobj %>% filter(!is.na(!!sym(case$group.by)))
        }, allow_empty)
        if (is.null(args$object)) { return() }
    }
    Idents(args$object) <- case$group.by

    cached <- get_cached(args, "FindAllMarkers", cache)
    if (!is.null(cached$data)) {
        log_info("  Using cached markers ...")
        markers <- cached$data
    } else {
        markers <- find_markers(args, find_all = TRUE)
        cached$data <- markers
        save_to_cache(cached, "FindAllMarkers", cache)
    }

    if (is.null(case$dotplot$assay)) {
        case$dotplot$assay <- case$assay
    }

    if (nrow(markers) == 0) {
        idents <- unique(Idents(args$object))
    } else {
        idents <- unique(markers$cluster)
    }
    for (ident in idents) {
        log_debug("  * Dealing with ident: {ident}...")
        if (case$prefix_group) {
            key <- paste0(info$section, "::", case$group.by, " - ", ident)
        } else {
            key <- paste0(info$section, "::", ident)
        }
        info_ident <- casename_info(key, cases, outdir, create = TRUE)
        if (nrow(markers) > 0) {
            markers_ident <- markers %>% filter(cluster == ident)
        } else {
            markers_ident <- markers
        }
        siggenes <- do_enrich(info_ident, markers_ident, case$sigmarkers, case$volcano_genes)

        if (length(siggenes) > 0) {
            args$ident.1 <- as.character(ident)
            do_dotplot(info_ident, siggenes, case$dotplot, args)
        }

        add_case_report(info_ident, case$sigmarkers, siggenes)
    }
}

find_markers <- function(findmarkers_args, find_all = FALSE) {
    if (find_all) {
        fun <- FindAllMarkers
        empty <- data.frame(
            gene = character(),
            p_val = numeric(),
            avg_log2FC = numeric(),
            pct.1 = numeric(),
            pct.2 = numeric(),
            p_val_adj = numeric(),
            cluster = character()
        )
    } else {
        fun <- FindMarkers
        empty <- data.frame(
            gene = character(),
            p_val = numeric(),
            avg_log2FC = numeric(),
            pct.1 = numeric(),
            pct.2 = numeric(),
            p_val_adj = numeric()
        )
    }

    call_findmarkers <- function(fn, args) {
        if (find_all) {
            do_call(fn, args)
        } else {
            do_call(fn, args) %>% rownames_to_column("gene")
        }
    }
    markers <- tryCatch({
        call_findmarkers(fun, findmarkers_args)
    }, error = function(e) {
        if (!grepl("PrepSCTFindMarkers", e$message) && defassay == "SCT") {
            log_warn(paste0("  ! ", e$message))
        }
        empty
    })

    if (nrow(markers) == 0 && defassay == "SCT") {
        log_warn("  ! No markers found from SCT assay, trying recorrect_umi = FALSE")
        findmarkers_args$recorrect_umi <- FALSE
        markers <- tryCatch({
            call_findmarkers(fun, findmarkers_args)
        }, error = function(e) {
            log_warn(paste0("  ! ", e$message))
            empty
        })
    }

    markers
}

sections <- c()
do_case <- function(casename) {
    if (isTRUE(cases[[casename]]$findall)) {
        log_info("- Dealing with case: {casename} (all idents) ...")
        do_case_findall(casename)
        return()
    }
    log_info("- Dealing with case: {casename} ...")

    info <- casename_info(casename, cases, outdir, create = TRUE)
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
    allow_empty = startsWith(case$group.by, "..")
    if (!is.null(case$subset)) {
        args$object <- ensure_sobj({
            srtobj %>% filter(!!parse_expr(case$subset) & !is.na(!!sym(case$group.by)))
        }, allow_empty)
        if (is.null(args$object)) { return() }
    } else {
        args$object <- ensure_sobj({
            srtobj %>% filter(!is.na(!!sym(case$group.by)))
        }, allow_empty)
        if (is.null(args$object)) { return() }
    }

    args$assay <- case$assay
    args$group.by <- case$group.by
    args$ident.1 <- case$ident.1
    args$ident.2 <- case$ident.2
    if (is.null(args$ident.2)) {
        args$ident.2 <- ".rest"
        args$object <- args$object %>% mutate(
            !!sym(args$group.by) := if_else(
                !!sym(args$group.by) == args$ident.1,
                args$ident.1,
                args$ident.2
            )
        )
    } else {
        args$object <- args$object %>%
            filter(!!sym(args$group.by) %in% c(args$ident.1, args$ident.2))
    }
    # args$logfc.threshold <- args$logfc.threshold %||% 0
    # args$min.cells.group <- args$min.cells.group %||% 1
    # args$min.cells.feature <- args$min.cells.feature %||% 1
    # args$min.pct <- args$min.pct %||% 0

    markers <- find_markers(args)
    siggenes <- do_enrich(info, markers, case$sigmarkers, case$volcano_genes)

    if (length(siggenes) > 0) {
        case$dotplot$assay <- case$dotplot$assay %||% args$assay
        do_dotplot(info, siggenes, case$dotplot, args)
    }

    sections <<- union(sections, info$section)
    if (info$section %in% overlapping_sections) {
        overlaps[[info$section]] <<- overlaps[[info$section]] %||% list()
        overlaps[[info$section]][[info$case]] <<- siggenes %||% character()
    }

    add_case_report(info, case$sigmarkers, siggenes)
}

do_overlap <- function(section) {
    log_info("- Dealing with overlapping: {section}...")

    ov_args <- overlapping[[section]]
    ov_dir <- file.path(outdir, "OVERLAPPING", section)
    dir.create(ov_dir, showWarnings = FALSE, recursive = TRUE)

    ov_cases <- overlaps[[section]]
    if (length(ov_cases) < 2) {
        stop(sprintf("  Not enough cases for overlap: %s", section))
    }

    if (is.list(ov_args$venn) && length(ov_cases) > 4) {
        stop(paste0("  Too many cases (", length(ov_cases)," > 4) for venn plot for section: ", section))
    }
    if (is.list(ov_args$venn)) {
        venn_plot <- file.path(ov_dir, "venn.png")
        venn_p <- ggVennDiagram(ov_cases, label_percent_digit = 1) +
            scale_fill_distiller(palette = "Reds", direction = 1) +
            scale_x_continuous(expand = expansion(mult = .2))
        ov_args$venn$devpars$file <- venn_plot
        do.call(png, ov_args$venn$devpars)
        print(venn_p)
        dev.off()
    }

    df_markers <- fromList(ov_cases)
    #  A  B  MARKERS
    #  1  0  G1
    #  1  0  G2
    #  0  1  G3
    #  0  1  G4
    #  1  1  G5
    df_markers$MARKERS = Reduce(union, ov_cases)
    df_markers = df_markers %>%
        group_by(across(-MARKERS)) %>%
        summarise(MARKERS = paste0(MARKERS, collapse = ","), .groups = "drop")

    write.table(
        df_markers,
        file.path(ov_dir, "markers.txt"),
        sep = "\t",
        row.names = FALSE,
        col.names = TRUE,
        quote = FALSE
    )

    if (is.list(ov_args$upset)) {
        upset_plot <- file.path(ov_dir, "upset.png")
        if (nrow(df_markers) == 0) {
            upset_p <- ggplot() +
                theme_void() +
                ggtitle("No overlapping markers found") +
                # center the title, and make it red
                theme(plot.title = element_text(hjust = 0.5, color = "red"))
            ov_args$upset$devpars <- list(
                res = 100, height = 42, width = 400
            )
        } else {
            upset_p <- upset(fromList(ov_cases))
        }
        ov_args$upset$devpars$file <- upset_plot
        do.call(png, ov_args$upset$devpars)
        print(upset_p)
        dev.off()
    }

    add_report(
        list(
            title = "Venn Diagram",
            ui = "flat",
            contents = list(
                list(
                    kind = "img",
                    src = file.path(ov_dir, "venn.png")
                )
            )
        ),
        list(
            title = "UpSet Plot",
            ui = "flat",
            contents = list(
                list(
                    kind = "img",
                    src = file.path(ov_dir, "upset.png")
                )
            )
        ),
        list(
            title = "Marker Table",
            ui = "flat",
            contents = list(
                list(
                    kind = "table",
                    data = list(nrows = 100),
                    src = file.path(ov_dir, "markers.txt")
                )
            )
        ),
        h1 = "Overlapping Markers",
        h2 = section,
        ui = "tabs"
    )
}

sapply(sort(names(cases)), do_case)

unhit_overlaps <- setdiff(overlapping_sections, names(overlaps))
if (length(unhit_overlaps) > 0) {
    log_warn(paste0("- No sections found for overlapping analysis: ", paste(unhit_overlaps, collapse = ", ")))
    log_warn("  Available sections: ", paste(sections, collapse = ", "))
}
sapply(sort(names(overlaps)), do_overlap)

save_report(joboutdir)
