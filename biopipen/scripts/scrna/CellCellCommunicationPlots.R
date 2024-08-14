{{ biopipen_dir | joinpaths: "utils", "misc.R" | source_r }}
library(rlang)
library(dplyr)
library(ggplot2)
library(CCPlotR)
{{ biopipen_dir | joinpaths: "scripts", "scrna", "CCPlotR-patch.R" | source_r }}

cccfile <- {{ in.cccfile | r }}
expfile <- {{ in.expfile | r }}
outdir <- {{ out.outdir | r }}
joboutdir <- {{ job.outdir | r }}
score_col <- {{ envs.score_col | r }}
subset <- {{ envs.subset | r }}
cases <- {{ envs.cases | r }}

ccc <- read.table(cccfile, header=TRUE, sep="\t", check.names = FALSE)
if (!is.null(subset)) {
    ccc <- ccc %>% dplyr::filter(!!parse_expr(subset))
}
if (ncol(ccc) > 10) {
    # from CellCellCommunication
    if (!is.null(expfile)) {
        log_warn("in.cccfile is from CellCellCommunication, in.expfile will be ignored")
    }
    if (is.null(score_col)) {
        stop("'envs.score_col' is required for CellCellCommunication output")
    }
    if (!score_col %in% colnames(ccc)) {
        stop(paste("Score column", score_col, "not found in the in.cccfile"))
    }
    # compose the expression data frame
    exp <- data.frame(
        cell_type = c(ccc$source, ccc$target),
        gene = c(ccc$ligand, ccc$receptor),
        mean_exp = c(ccc$ligand_trimean, ccc$receptor_trimean)
    ) %>% distinct()
    ccc <- ccc %>% select(
        source, target,
        ligand, receptor,
        !!sym(score_col)
    ) %>% rename(score = !!sym(score_col))
} else {
    if (!is.null(expfile)) {
        exp <- read.table(expfile, header=TRUE, sep="\t", check.names = FALSE)
    }
}

if (length(cases) == 0) {
    stop("No cases provided.")
}

.get_default_devpars <- function(kind, nrows, ncols = NULL) {
    if (kind == "arrow") {
        list(
            res = 100,
            width = 600,
            height = 50 + nrows * 20
        )
    } else if (kind == "circos") {
        list(
            res = 100,
            width = 800,
            height = 800
        )
    } else if (kind == "dotplot") {
        list(
            res = 100,
            width = 120 + ncols * 60,
            height = 300 + nrows * 40
        )
    } else if (kind == "heatmap") {
        list(
            res = 100,
            width = 120 + ncols * 60,
            height = 300 + ncols * 40
        )
    } else if (kind == "network") {
        list(
            res = 100,
            width = 1200,
            height = 1200
        )
    } else if (kind == "sigmoid") {
        list(
            res = 100,
            width = max(800, ncols * 200),
            height = 100 + nrows * 60
        )
    }
}

images <- lapply(names(cases), function(name) {
    log_info("- Case: ", name, " ...")
    case <- cases[[name]]

    kind <- match.arg(case$kind, c("arrow", "circos", "dotplot", "heatmap", "network", "sigmoid"))
    fun <- get(paste0("cc_", kind))
    case$kind <- NULL

    gg <- NULL
    if (kind == "arrow") {
        cell_types <- case$cell_types
        if (is.null(cell_types) || length(cell_types) != 2) {
            stop("'case.cell_types' is required and must be a vector of length 2")
        }
        n_ligand <- length(unique(ccc[ccc$source == cell_types[1], "ligand"]))
        n_receptor <- length(unique(ccc[ccc$target == cell_types[2], "receptor"]))
        default_devpars <- .get_default_devpars(kind, nrows = max(n_ligand, n_receptor))
    } else if (kind == "circos") {
        nrows <- length(unique(c(ccc$source, ccc$target)))
        default_devpars <- .get_default_devpars(kind, nrows = nrows)
    } else if (kind == "dotplot" || kind == "heatmap") {
        nrows <- length(unique(ccc$source))
        ncols <- length(unique(ccc$target))
        default_devpars <- .get_default_devpars(kind, nrows = nrows, ncols = ncols)
        if (
            (kind == "heatmap" && (is.null(case$option) || case$option != "B")) ||
            (kind == "dotplot" && (is.null(case$option) || case$option != "B"))) {
            gg <- theme(axis.text.x = element_text(angle = 90, hjust = 1))
        }
    } else if (kind == "network") {
        nrows <- length(unique(c(ccc$source, ccc$target)))
        ncols <- length(unique(c(ccc$ligand, ccc$receptor)))
        default_devpars <- .get_default_devpars(kind, nrows = nrows, ncols = ncols)
        gg <- theme(plot.margin = margin(c(50, 50, 50, 50), "pt"))
    } else if (kind == "sigmoid") {
        nrows <- (case$n_top_ints %||% 20) / 2  # approx
        ncols <- length(unique(c(ccc$source, ccc$target))) / 2
        default_devpars <- .get_default_devpars(kind, nrows = nrows, ncols = ncols)
    }
    devpars <- case$devpars %||% default_devpars
    devpars$res <- devpars$res %||% default_devpars$res
    devpars$width <- devpars$width %||% default_devpars$width
    devpars$height <- devpars$height %||% default_devpars$height
    case$devpars <- NULL

    section <- case$section
    case$section <- NULL

    case$cc_df <- ccc
    if ("exp_df" %in% names(formals(fun))) {
        case$exp_df <- exp
    }
    outpath <- file.path(outdir, paste0(slugify(name), ".png"))
    png(outpath, width=devpars$width, height=devpars$height, res=devpars$res)
    p <- do_call(fun, case)
    if (!is.null(gg)) { p <- p + gg }
    print(p)
    dev.off()

    list(
        section = section,
        kind = "table_image",
        src = outpath,
        name = name
    )
})

section_images = list()
for (image in images) {
    section <- image$section
    image$section <- NULL
    if (is.null(section)) {
        section = "DEFAULT"
    }
    if (!section %in% names(section_images)) {
        section_images[[section]] = list()
    }
    section_images[[section]][[length(section_images[[section]]) + 1]] = image
}

if (length(section_images) == 1 && names(section_images)[1] == "DEFAULT") {
    add_report(
        section_images,
        h1 = "Cell-Cell Communication Plots",
        ui = "table_of_images"
    )
} else {
    for (section in names(section_images)) {
        imgplots = section_images[[section]]
        add_report(
            list(
                ui = "table_of_images",
                contents = imgplots
            ),
            h1 = section
        )
    }
}

save_report(joboutdir)
