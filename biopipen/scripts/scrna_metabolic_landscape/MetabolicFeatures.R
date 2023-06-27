source("{{biopipen_dir}}/utils/gsea.R")

library(parallel)

sobjfile <- {{ in.sobjfile | r }}
outdir <- {{ out.outdir | r }}
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

pathways <- gmt_pathways(gmtfile)
metabolics <- unique(as.vector(unname(unlist(pathways))))
sobj <- readRDS(sobjfile)

do_one_group <- function(obj, group, outputdir) {
    print(paste("- Processing group", grouping, ":", group))
    groupname = paste0(grouping_prefix, group)
    odir = file.path(outputdir, groupname)
    dir.create(odir, showWarnings = FALSE)

    classes = as.character(obj@meta.data[[grouping]])
    classes[classes != group] <- "_REST"
    classes[classes == group] <- groupname
    exprs = GetAssayData(obj)
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
    }, error=function(e) {
        unlink(odir, recursive = T, force = T)
        warning(paste("Unable to run for", group))
        warning(e)
    })

}

do_one_subset <- function(s, subset_col, subset_prefix) {
    print(paste("Processing subset", subset_col, ":", s))
    if (is.null(s)) {
        outputdir <- file.path(outdir, "ALL")
        subset_obj <- sobj
    } else {
        outputdir <- file.path(outdir, paste0(subset_prefix, s))
        subset_code <- paste0("subset(sobj, subset = ", subset_col, "=='", s, "')")
        subset_obj <- eval(parse(text = subset_code))
    }
    dir.create(outputdir, showWarnings = FALSE)

    subset_obj <- subset(subset_obj, features = intersect(rownames(subset_obj), metabolics))

    groups = subset_obj@meta.data[[grouping]]
    x = mclapply(as.character(unique(groups)), function(group) {
        do_one_group(subset_obj, group, outputdir)
    }, mc.cores = ncores)
    if (any(unlist(lapply(x, class)) == "try-error")) {
        stop("mclapply error")
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

