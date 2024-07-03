source("{{biopipen_dir}}/utils/misc.R")
source("{{biopipen_dir}}/utils/gsea.R")

library(parallel)
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

do_one_group <- function(obj, features, group, outputdir, h1) {
    log_info(paste("- Processing group", grouping, ":", group))
    groupname = paste0(grouping_prefix, group)
    odir = file.path(outputdir, slugify(groupname))
    dir.create(odir, showWarnings = FALSE)

    classes = as.character(obj@meta.data[[grouping]])
    classes[classes != group] <- "_REST"
    classes[classes == group] <- groupname
    if (any(table(classes) < 5)) {
        msg <- paste("  Skipped. One of the groups has less than 5 cells.")
        log_warn(msg)
        # write a warning.txt to odir with the message and table(classes)
        write(paste0(msg, "\n\n"), file = file.path(odir, "warning.txt"))
        write.table(
            table(classes),
            file = file.path(odir, "warning.txt"),
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            append = TRUE
        )
        return(
            list(
                list(kind = "error", content = msg),
                h1 = ifelse(is.null(h1), groupname, h1),
                h2 = ifelse(is.null(h1), "#", groupname)
            )
        )
    }

    exprs = GetAssayData(obj)[features, , drop = FALSE]
    tryCatch({
        if (fgsea) {
            runFGSEA(
                prerank(exprs, groupname, "_REST", classes, method = prerank_method),
                gmtfile,
                top = top,
                outdir = odir
            )
        } else {
            runGSEA(
                exprs,
                classes,
                gmtfile,
                odir
            )
        }

        # Can't add report directly, mclapply can't modify global variables.
        report = list(
            list(kind = "fgsea", dir = odir),
            h1 = ifelse(is.null(h1), groupname, h1),
            h2 = ifelse(is.null(h1), "#", groupname)
        )
    }, error=function(e) {
        unlink(odir, recursive = T, force = T)
        log_warn(paste("Unable to run for", group))
        log_warn(e$message)

        report = list(
            list(
                kind = "error",
                content = paste0("Error running GSEA for ", group, ": ", e$message)
            ),
            h1 = ifelse(is.null(h1), groupname, h1),
            h2 = ifelse(is.null(h1), "#", groupname)
        )
    })

    report
}

do_one_subset <- function(s, subset_col, subset_prefix) {
    log_info(paste("Processing subset", subset_col, ":", s))
    if (is.null(s)) {
        outputdir <- file.path(outdir, "ALL")
        subset_obj <- sobj
    } else {
        outputdir <- file.path(outdir, slugify(paste0(subset_prefix, s)))
        subset_code <- paste0("subset(sobj, subset = ", subset_col, "=='", s, "')")
        subset_obj <- eval(parse(text = subset_code))
    }
    dir.create(outputdir, showWarnings = FALSE)

    # subset_obj <- subset(subset_obj, features = intersect(rownames(subset_obj), metabolics))
    features = intersect(rownames(subset_obj), metabolics)

    h1 <- NULL
    if (!is.null(s)) {
        h1 <- paste0(subset_prefix, s)
    }
    groups = subset_obj@meta.data[[grouping]]
    x = mclapply(as.character(unique(groups)), function(group) {
        do_one_group(subset_obj, features, group, outputdir, h1)
    }, mc.cores = ncores)
    if (any(unlist(lapply(x, class)) == "try-error")) {
        stop("mclapply error")
    }
    for (r in x) {
        if (!is.null(r)) {
            do.call(add_report, r)
        }
    }
}

do_one_subset_col <- function(subset_col, subset_prefix) {
    if (is.null(subset_col)) {
        do_one_subset(NULL, subset_col = NULL, subset_prefix = NULL)
    }
    subsets <- na.omit(unique(sobj@meta.data[[subset_col]]))
    lapply(subsets, do_one_subset, subset_col = subset_col, subset_prefix = subset_prefix)
}

if (is.null(subsetting_cols)) {
    do_one_subset_col(NULL)
} else {
    for (i in seq_along(subsetting_cols)) {
        do_one_subset_col(subsetting_cols[i], subsetting_prefix[i])
    }
}

save_report(joboutdir)
