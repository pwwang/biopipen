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

log_info("Setting up EnrichR ...")
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
assay <- {{ envs.assay | r }}
sigmarkers <- {{ envs.sigmarkers | r }}
volcano_genes <- {{ envs.volcano_genes | r }}
subset <- {{ envs.subset | r }}
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
    log_info("- Running enrichment for case: {case}")
    parts <- strsplit(case, ":")[[1]]
    sec <- parts[1]
    case <- paste0(parts[-1], collapse = ":")
    casedir <- file.path(outdir, sec, case)
    dir.create(casedir, showWarnings = FALSE, recursive = TRUE)
    if (nrow(markers) == 0) {
        log_warn("  No markers found for case: {case}")
        cat("No markers found.", file = file.path(casedir, "error.txt"))
        return()
    }
    plot_volcano(markers, file.path(casedir, "volcano.png"), sig, volgenes)
    markers_sig <- markers %>% filter(!!parse_expr(sig))
    if (nrow(markers_sig) == 0) {
        log_warn("  No significant markers found for case: {case}")
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
    log_info("Dealing with case: {casename}...")
    sec_case_names <- strsplit(casename, ":")[[1]]
    cname <- paste(sec_case_names[-1], collapse = ":")
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
        args$object <- srtobj %>% filter(!!parse_expr(case$subset) & filter(!is.na(!!sym(case$group.by))))
    } else {
        args$object <- srtobj %>% filter(!is.na(!!sym(case$group.by)))
    }
    markers <- tryCatch({
        do_call(FindMarkers, args) %>% rownames_to_column("gene")
    }, error = function(e) {
        warning(e$message, immediate. = TRUE)
        data.frame(
            gene = character(),
            p_val = numeric(),
            avg_log2FC = numeric(),
            pct.1 = numeric(),
            pct.2 = numeric(),
            p_val_adj=numeric()
        )
    })
    do_enrich(casename, markers, case$sigmarkers, case$volcano_genes)

    siggenes <- markers %>%
        filter(!!parse_expr(case$sigmarkers)) %>%
        pull(gene) %>%
        unique()

    if (length(siggenes) > 0) {
        dotplot_devpars <- case$dotplot$devpars
        if (is.null(args$ident.2)) {
            case$dotplot$object <- args$object
            case$dotplot$object@meta.data <- case$dotplot$object@meta.data %>%
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
            case$dotplot$object <- args$object %>%
                filter(!!sym(args$group.by) %in% c(args$ident.1, args$ident.2)) %>%
                mutate(!!sym(args$group.by) := factor(
                    !!sym(args$group.by),
                    levels = c(args$ident.1, args$ident.2)
                ))
        }
        case$dotplot$devpars <- NULL
        case$dotplot$features <- siggenes
        case$dotplot$group.by <- args$group.by
        case$dotplot$assay <- case$assay
        dotplot_width = ifelse(
            is.null(dotplot_devpars$width),
            if (length(siggenes) <= 20) length(siggenes) * 60 else length(siggenes) * 30,
            dotplot_devpars$width
        )
        dotplot_height = ifelse(is.null(dotplot_devpars$height), 600, dotplot_devpars$height)
        dotplot_res = ifelse(is.null(dotplot_devpars$res), 100, dotplot_devpars$res)
        dotplot_file <- file.path(outdir, sec_case_names[1], cname, "dotplot.png")
        png(dotplot_file, res = dotplot_res, width = dotplot_height, height = dotplot_width)
        # rotate x axis labels
        print(
            do_call(DotPlot, case$dotplot) +
            theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
            coord_flip()
        )
        dev.off()
    }

    if (sec_case_names[1] %in% overlap) {
        if (is.null(overlaps[[sec_case_names[1]]])) {
            overlaps[[sec_case_names[1]]] <<- list()
        }
        overlaps[[sec_case_names[1]]][[cname]] <<- siggenes
    }
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
}

sapply(sort(names(cases)), do_case)
sapply(sort(names(overlaps)), do_overlap)
