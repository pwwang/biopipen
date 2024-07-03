source("{{biopipen_dir}}/utils/misc.R")
source("{{biopipen_dir}}/utils/gsea.R")

library(parallel)
library(scater)
library(Seurat)

sobjfile <- {{ in.sobjfile | r }}
outdir <- {{ out.outdir | r }}
joboutdir <- {{ job.outdir | r }}
gmtfile <- {{ envs.gmtfile | r }}
ncores <- {{ envs.ncores | r }}
fgsea <- {{ envs.fgsea | r }}
top <- {{ envs.top | r }}
prerank_method <- {{ envs.prerank_method | r }}
grouping <- {{ envs.grouping | r }}
grouping_prefix <- {{ envs.grouping_prefix | r }}
subsetting_cols <- {{ envs.subsetting | r }}
subsetting_prefix <- {{ envs.subsetting_prefix | r }}
subsetting_comparison <- {{ envs.subsetting_comparison | r }}

if (!is.null(grouping_prefix) && nchar(grouping_prefix) > 0) {
    grouping_prefix = paste0(grouping_prefix, "_")
}

if (!is.null(subsetting_prefix) && nchar(subsetting_prefix) > 0) {
    subsetting_prefix = paste0(subsetting_prefix, "_")
}

set.seed(8525)

## gmt_pathways is copied from fgsea package.
gmt_pathways <- function(gmt_file) {
    pathway_lines <- strsplit(readLines(gmt_file), "\t")
    pathways <- lapply(pathway_lines, tail, -2)
    names(pathways) <- sapply(pathway_lines, head, 1)
    pathways
}

gmtfile <- localizeGmtfile(gmtfile)
pathways <- gmt_pathways(gmtfile)
metabolics <- unique(as.vector(unname(unlist(pathways))))
sobj <- readRDS(sobjfile)

do_one_comparison <- function(
    obj,
    compname,
    genes,
    case,
    control,
    groupdir,
    subset_col,
    subset_prefix,
    groupname
) {
    log_info(paste("  Design: {compname} ({case}, {control})"))
    case_code = paste0("subset(obj, subset = ", subset_col, " == '", case, "')")
    case_obj = tryCatch({
        eval(parse(text = case_code))
    }, error = function(e) {
        NULL
    })
    if (is.null(case_obj)) {
        log_warn("          Skip (not enough cells in case)")
        return (NULL)
    }
    control_code = paste0("subset(obj, subset = ", subset_col, " == '", control, "')")
    control_obj = tryCatch({
        eval(parse(text = control_code))
    }, error = function(e) {
        NULL
    })
    if (is.null(control_obj)) {
        log_warn("          Skip (not enough cells in control)")
        add_report(
            list(kind = "error", content = "Not enough cells in control"),
            h1 = groupname,
            h2 = compname
        )
        return (NULL)
    }
    exprs_case = GetAssayData(case_obj)[genes, , drop = FALSE]
    exprs_control = GetAssayData(control_obj)[genes, , drop = FALSE]

    odir = file.path(groupdir, paste0(subset_prefix, compname))
    dir.create(odir, showWarnings = FALSE)
    if (ncol(exprs_case) < 5 || ncol(exprs_control) < 5) {
        log_warn("  Skipped (not enough cells).")
        wfile <- file.path(odir, "warning.txt")
        write("Skipped (not enough cells)\n\n", file = wfile)
        write(paste0("n_cells (Case):", ncol(exprs_case)), file = wfile, append = TRUE)
        write(paste0("n_cells (Control):", ncol(exprs_control)), file = wfile, append = TRUE)

        return(list(
            list(kind = "error", content = "Not enough cells"),
            h1 = groupname,
            h2 = compname
        ))
    }
    if (fgsea) {
        ranks = prerank(
            cbind(exprs_case, exprs_control),
            case,
            control,
            c(rep(case, ncol(exprs_case)), rep(control, ncol(exprs_control))),
            method = prerank_method
        )

        runFGSEA(
            ranks,
            gmtfile,
            top = top,
            outdir = odir,
            envs = list(nproc = 1)
        )

        report = list(
            list(kind = "fgsea", dir = odir),
            h1 = groupname,
            h2 = compname
        )
    } else {
        runGSEA(
            cbind(exprs_case, exprs_control),
            c(rep(case, ncol(exprs_case)), rep(control, ncol(exprs_control))),
            gmtfile,
            odir
        )

        report = list()
    }

    report
}

do_one_group <- function(group) {
    log_info("- Group: {group} ...")

    genes = intersect(metabolics, rownames(sobj))
    group_code = paste0(
        "subset(sobj, subset = ", grouping, " == '", group, "')"
    )
    obj = eval(parse(text = group_code))
    groupname = paste0(grouping_prefix, group)
    groupdir = file.path(outdir, slugify(groupname))
    dir.create(groupdir, showWarnings = FALSE)

    report = list()
    for (i in seq_along(subsetting_comparison)) {
        sci = subsetting_comparison[[i]]
        if (is.null(sci) || length(sci) == 0) {
            next
        }
        rs = lapply(
            names(sci),
            function(compname) {
                do_one_comparison(
                    obj,
                    compname,
                    genes,
                    sci[[compname]][1],
                    sci[[compname]][2],
                    groupdir,
                    subsetting_cols[i],
                    subsetting_prefix[i],
                    groupname
                )
            }
        )
        if (length(rs) > 0) {
            report = c(report, rs)
        }
    }
    report
}

groups = sort(as.character(unique(sobj@meta.data[[grouping]])))
if (ncores == 1) {
    x = lapply(groups, do_one_group)
} else {
    x = mclapply(groups, do_one_group, mc.cores = ncores)
    if (any(unlist(lapply(x, class)) == "try-error")) {
        stop("mclapply error")
    }
}
report = unlist(x, recursive = FALSE)
for (r in report) {
    if (!is.null(r)) {
        do.call(add_report, r)
    }
}

save_report(joboutdir)
