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
library(ggVennDiagram)
library(UpSetR)
library(slugify)

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
section <- {{ envs.section | r }}
dbs <- {{ envs.dbs | r }}
assay <- {{ envs.assay | r }}
sigmarkers <- {{ envs.sigmarkers | r }}
volcano_genes <- {{ envs.volcano_genes | r }}
subset <- {{ envs.subset | r }}
use_presto <- {{ envs.use_presto | r }}
rest <- {{ envs.rest | r: todot="-" }}
dotplot <- {{ envs.dotplot | r: todot="-" }}
cases <- {{ envs.cases | r: todot="-" }}
overlap <- {{ envs.overlap | r }}

overlaps <- list()

if (is.character(volcano_genes) && length(volcano_genes) == 1) {
    volcano_genes <- trimws(strsplit(volcano_genes, ",")[[1]])
}

set.seed(8525)
if (ncores > 1) {
    options(future.globals.maxSize = 80000 * 1024^2)
    plan(strategy = "multicore", workers = ncores)
}

log_info("Reading Seurat object ...")
srtobj <- readRDS(srtfile)

log_info("Mutate meta data if needed ...")
if (!is.null(mutaters) && length(mutaters) > 0) {
    srtobj@meta.data <- srtobj@meta.data %>%
        mutate(!!!lapply(mutaters, parse_expr))
}

log_info("Expanding cases ...")
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
            assay = assay,
            subset = subset,
            use_presto = use_presto,
            sigmarkers = sigmarkers,
            volcano_genes = volcano_genes,
            dotplot = dotplot,
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
            assay = assay,
            subset = subset,
            use_presto = use_presto,
            sigmarkers = sigmarkers,
            volcano_genes = volcano_genes,
            dotplot = dotplot,
            rest = rest
        )
        case$rest <- list_update(rest, case$rest)
        case$dotplot$devpars <- list_update(dotplot$devpars, case$dotplot$devpars)
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

sections <- c()

newcases <- list()
for (name in names(cases)) {
    case <- cases[[name]]
    if (is.null(case$each) && !is.null(case$ident.1)) {
        sections <- c(sections, case$section)
        newcases[[paste0(case$section, ":", name)]] <- case
    } else if (is.null(case$each)) {
        # is.null(case$ident.1)
        sections <- c(sections, name)
        newcases[[name]] <- case
        newcases[[name]]$findall <- TRUE
        # idents <- srtobj@meta.data %>% pull(case$group.by) %>% unique() %>% na.omit()
        # for (ident in idents) {
        #     newcases[[paste0(name, ":", ident)]] <- case
        #     newcases[[paste0(name, ":", ident)]]$ident.1 <- ident
        # }
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
                kname <- if (name == "DEFAULT") "" else paste0(" - ", name)
                sections <- c(sections, paste0(each, kname))
                key <- paste0(each, kname)
                newcases[[key]] <- case
                newcases[[key]]$group.by <- by
                newcases[[key]]$findall <- TRUE
                # idents <- srtobj@meta.data %>% pull(case$group.by) %>% unique() %>% na.omit()
                # for (ident in idents) {
                #     key <- paste0(each, kname, ":", ident)
                #     if (case$prefix_each) {
                #         key <- paste0(case$each, " - ", key)
                #     }
                #     newcases[[key]] <- case
                #     newcases[[key]]$ident.1 <- ident
                #     newcases[[key]]$group.by <- by
                # }
            } else {
                sections <- c(sections, case$each)
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
single_section <- length(unique(sections)) == 1

casename_info <- function(casename, create = FALSE) {
    sec_case_names <- strsplit(casename, ":")[[1]]
    cname <- paste(sec_case_names[-1], collapse = ":")

    out <- list(
        casename = casename,
        section = sec_case_names[1],
        case = cname,
        section_slug = slugify(sec_case_names[1], tolower = FALSE),
        case_slug = slugify(cname, tolower = FALSE)
    )
    out$casedir <- file.path(outdir, out$section_slug, out$case_slug)
    if (create) {
        dir.create(out$casedir, showWarnings = FALSE, recursive = TRUE)
    }
    out
}

plot_volcano = function(markers, volfile, sig, volgenes) {
    # markers
    #                  gene        p_val avg_log2FC pct.1 pct.2    p_val_adj
    # 1            CCL5 1.883596e-11 -4.8282535 0.359 0.927 4.332270e-09
    # 2        HLA-DQB1 3.667713e-09  6.1543174 0.718 0.098 8.435740e-07
    # 3        HLA-DRB5 1.242993e-07  3.9032231 0.744 0.195 2.858885e-05
    # 4           CD79B 2.036731e-07  4.2748835 0.692 0.146 4.684482e-05
    log_info("- Plotting volcano plot ...")
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
    log_info("- Running enrichment for case: {info$casename}")

    if (nrow(markers) == 0) {
        log_warn("  No markers found for case: {info$casename}")
        return(NULL)
    }

    plot_volcano(markers, file.path(info$casedir, "volcano.png"), sig, volgenes)

    markers_sig <- markers %>% filter(!!parse_expr(sig))
    if (nrow(markers_sig) == 0) {
        log_warn("  No significant markers found for case: {info$casename}")
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
        return(NULL)
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
    dotplot_devpars <- dotplot$devpars
    if (is.null(args$ident.2)) {
        dotplot$object <- args$object
        dotplot$object@meta.data <- dotplot$object@meta.data %>%
            mutate(
                !!sym(args$group.by) := if_else(
                    !!sym(args$group.by) == args$ident.1,
                    args$ident.1,
                    ".Other"
                ),
                !!sym(args$group.by) := factor(
                    !!sym(args$group.by),
                    levels = c(args$ident.1, ".Other")
                )
            )
    } else {
        dotplot$object <- args$object %>%
            filter(!!sym(args$group.by) %in% c(args$ident.1, args$ident.2)) %>%
            mutate(!!sym(args$group.by) := factor(
                !!sym(args$group.by),
                levels = c(args$ident.1, args$ident.2)
            ))
    }
    dotplot$devpars <- NULL
    dotplot$features <- siggenes
    dotplot$group.by <- args$group.by
    dotplot_width = ifelse(
        is.null(dotplot_devpars$width),
        if (length(siggenes) <= 20) length(siggenes) * 60 else min(1000, length(siggenes)) * 30,
        dotplot_devpars$width
    )
    dotplot_height = ifelse(is.null(dotplot_devpars$height), 600, dotplot_devpars$height)
    dotplot_res = ifelse(is.null(dotplot_devpars$res), 100, dotplot_devpars$res)
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
    h1 = ifelse(
        info$section == "DEFAULT",
        info$case,
        ifelse(
            single_section,
            paste0(
                ifelse(info$section == "seurat_clusters", "Cluster", info$section),
                " - ",
                info$case
            ),
            info$section
        )
    )
    h2 = ifelse(
        info$section == "DEFAULT",
        "#",
        ifelse(single_section, "#", info$case)
    )
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
    if (is.null(siggenes)) {
        add_report(
            list(
                kind = "error",
                content = "No enough significant markers found for enrichment analysis"
            ),
            h1 = h1,
            h2 = ifelse(h2 == "#", "Enrichment Analysis", h2),
            h3 = ifelse(h2 == "#", "#", "Enrichment Analysis"),
            ui = "flat"
        )
    } else {
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


do_case_findall <- function(casename) {
    log_info("- Using FindAllMarkers for case: {casename}...")

    case <- cases[[casename]]
    if (case$use_presto) {
        log_info("- Using presto::wilcoxauc for case: {casename}...")
        args = list(
            object = if (!is.null(case$subset)) {
                srtobj %>% filter(!!parse_expr(case$subset) & !is.na(!!sym(case$group.by)))
            } else {
                srtobj %>% filter(!is.na(!!sym(case$group.by)))
            },
            group.by = case$group.by,
            ident.1 = case$ident.1,
            ident.2 = case$ident.2
        )
        markers <- tryCatch({
            # FindAllMarkers:
            # gene, p_val, avg_log2FC, pct.1, pct.2, p_val_adj, cluster
            # wilcoxauc:
            # feature, group, avgExpr, logFC, statistic, auc, pval, padj, pct_in, pct_out
            presto::wilcoxauc(
                args$object,
                group_by = case$group.by,
                seurat_assay = case$assay
            ) %>% select(
                gene = feature,
                p_val = pval,
                avg_log2FC = logFC,
                pct.1 = pct_in,
                pct.2 = pct_out,
                p_val_adj = padj,
                cluster = group
            )
        }, error = function(e) {
            log_warn(e$message)
            data.frame(
                gene = character(),
                p_val = numeric(),
                avg_log2FC = numeric(),
                pct.1 = numeric(),
                pct.2 = numeric(),
                p_val_adj=numeric(),
                cluster = character()
            )
        })
    } else {
        args <- case$rest
        args$group.by <- case$group.by
        if (is.null(args$logfc.threshold)) {
            args$locfc.threshold <- 0
        }
        if (is.null(args$min.cells.group)) {
            args$min.cells.group <- 1
        }
        if (is.null(args$min.cells.feature)) {
            args$min.cells.feature <- 1
        }
        if (is.null(args$min.pct)) {
            args$min.pct <- 0
        }
        if (!is.null(case$subset)) {
            args$object <- srtobj %>% filter(!!parse_expr(case$subset) & !is.na(!!sym(case$group.by)))
        } else {
            args$object <- srtobj %>% filter(!is.na(!!sym(case$group.by)))
        }
        Idents(args$object) <- case$group.by

        markers <- tryCatch({
            do_call(FindAllMarkers, args)
            # gene, p_val, avg_log2FC, pct.1, pct.2, p_val_adj, cluster
        }, error = function(e) {
            log_warn(e$message)

            data.frame(
                gene = character(),
                p_val = numeric(),
                avg_log2FC = numeric(),
                pct.1 = numeric(),
                pct.2 = numeric(),
                p_val_adj=numeric(),
                cluster = character()
            )
        })

        if (nrow(markers) == 0 && DefaultAssay(srtobj) == "SCT") {
            log_warn("  No markers found from SCT assay, try recorrect_umi = FALSE")
            args$recorrect_umi <- FALSE
            markers <- tryCatch({
                do_call(FindAllMarkers, args)
            }, error = function(e) {
                log_warn(e$message)
                data.frame(
                    gene = character(),
                    p_val = numeric(),
                    avg_log2FC = numeric(),
                    pct.1 = numeric(),
                    pct.2 = numeric(),
                    p_val_adj=numeric(),
                    cluster = character()
                )
            })
        }
    }

    if (is.null(case$dotplot$assay)) {
        case$dotplot$assay <- assay
    }
    idents <- unique(markers$cluster)
    for (ident in idents) {
        log_info("- Dealing with ident: {ident}...")
        info <- casename_info(paste0(casename, ":", ident), create = TRUE)
        siggenes <- do_enrich(info, markers %>% filter(cluster == ident), case$sigmarkers, case$volcano_genes)

        if (length(siggenes) > 0) {
            args$ident.1 <- as.character(ident)
            do_dotplot(info, siggenes, case$dotplot, args)
        }
        add_case_report(info, case$sigmarkers, siggenes)

        if (info$section %in% overlap) {
            if (is.null(overlaps[[info$section]])) {
                overlaps[[info$section]] <<- list()
            }
            overlaps[[info$section]][[info$case]] <<- siggenes
        }
    }
}


do_case <- function(casename) {
    log_info("Dealing with case: {casename}...")

    if (isTRUE(cases[[casename]]$findall)) {
        do_case_findall(casename)
        return()
    }

    info <- casename_info(casename, create = TRUE)
    case <- cases[[casename]]
    # ident1
    # ident2
    # groupby
    # each  # expanded
    # prefix_each
    # dbs
    # sigmarkers
    # rest
    if (case$use_presto) {
        log_info("- Using presto::wilcoxauc for case: {casename}...")
        args = list()
        markers <- tryCatch({
            if (!is.null(case$subset)) {
                args$object <- srtobj %>% filter(!!parse_expr(case$subset) & !is.na(!!sym(case$group.by)))
            } else {
                args$object <- srtobj %>% filter(!is.na(!!sym(case$group.by)))
            }
            if (!is.null(case$ident.2)) {
                args$object <- args$object %>% filter(!!sym(case$group.by) %in% c(case$ident.1, case$ident.2))
            }
            args$group.by <- case$group.by
            args$ident.1 <- case$ident.1
            args$ident.2 <- case$ident.2
            # FindMarkers:
            # gene, p_val, avg_log2FC, pct.1, pct.2, p_val_adj
            # wilcoxauc:
            # feature, group, avgExpr, logFC, statistic, auc, pval, padj, pct_in, pct_out
            presto::wilcoxauc(
                args$object,
                group_by = case$group.by,
                seurat_assay = case$assay
            ) %>%
            filter(group == case$ident.1) %>%
            select(
                gene = feature,
                p_val = pval,
                avg_log2FC = logFC,
                pct.1 = pct_in,
                pct.2 = pct_out,
                p_val_adj = padj
            )
        }, error = function(e) {
            log_warn(e$message)
            data.frame(
                gene = character(),
                p_val = numeric(),
                avg_log2FC = numeric(),
                pct.1 = numeric(),
                pct.2 = numeric(),
                p_val_adj=numeric()
            )
        })
    } else {
        args <- case$rest
        args$group.by <- case$group.by
        args$ident.1 <- case$ident.1
        args$ident.2 <- case$ident.2
        if (is.null(args$logfc.threshold)) {
            args$locfc.threshold <- 0
        }
        if (is.null(args$min.cells.group)) {
            args$min.cells.group <- 1
        }
        if (is.null(args$min.cells.feature)) {
            args$min.cells.feature <- 1
        }
        if (is.null(args$min.pct)) {
            args$min.pct <- 0
        }
        if (!is.null(case$subset)) {
            args$object <- srtobj %>% filter(!!parse_expr(case$subset) & !is.na(!!sym(case$group.by)))
        } else {
            args$object <- srtobj %>% filter(!is.na(!!sym(case$group.by)))
        }

        markers <- tryCatch({
            do_call(FindMarkers, args) %>% rownames_to_column("gene")
        }, error = function(e) {
            log_warn(e$message)
            data.frame(
                gene = character(),
                p_val = numeric(),
                avg_log2FC = numeric(),
                pct.1 = numeric(),
                pct.2 = numeric(),
                p_val_adj = numeric()
            )
        })

        if (nrow(markers) == 0 && DefaultAssay(srtobj) == "SCT") {
            log_warn("  No markers found from SCT assay, try recorrect_umi = FALSE")
            args$recorrect_umi <- FALSE
            markers <- tryCatch({
                do_call(FindMarkers, args) %>% rownames_to_column("gene")
            }, error = function(e) {
                log_warn(e$message)
                data.frame(
                    gene = character(),
                    p_val = numeric(),
                    avg_log2FC = numeric(),
                    pct.1 = numeric(),
                    pct.2 = numeric(),
                    p_val_adj=numeric(),
                    cluster = character()
                )
            })
        }
    }

    siggenes <- do_enrich(info, markers, case$sigmarkers, case$volcano_genes)

    if (length(siggenes) > 0) {
        if (is.null(case$dotplot$assay)) {
            case$dotplot$assay <- assay
        }
        do_dotplot(info, siggenes, case$dotplot, args)
    }

    if (info$section %in% overlap) {
        if (is.null(overlaps[[info$section]])) {
            overlaps[[info$section]] <<- list()
        }
        overlaps[[info$section]][[info$case]] <<- siggenes
    }

    add_case_report(info, case$sigmarkers, siggenes)
}

do_overlap <- function(section) {
    log_info("Dealing with overlap: {section}...")

    ov_dir <- file.path(outdir, "OVERLAPS", section)
    dir.create(ov_dir, showWarnings = FALSE, recursive = TRUE)

    ov_cases <- overlaps[[section]]
    if (length(ov_cases) < 2) {
        stop(sprintf("  Not enough cases for overlap: %s", section))
    }

    if (length(ov_cases) <= 4) {
        venn_plot <- file.path(ov_dir, "venn.png")
        venn_p <- ggVennDiagram(ov_cases, label_percent_digit = 1) +
            scale_fill_distiller(palette = "Reds", direction = 1) +
            scale_x_continuous(expand = expansion(mult = .2))
        png(venn_plot, res = 100, width = 1000, height = 600)
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

    upset_plot <- file.path(ov_dir, "upset.png")
    upset_p <- upset(fromList(ov_cases))
    png(upset_plot, res = 100, width = 800, height = 600)
    print(upset_p)
    dev.off()

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
sapply(sort(names(overlaps)), do_overlap)

save_report(joboutdir)
