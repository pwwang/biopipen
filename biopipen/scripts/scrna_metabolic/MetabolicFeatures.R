source("{{biopipen_dir}}/utils/gsea.R")

library(parallel)

sceobjfile <- {{ in.sceobj | r }}
gmtfile <- {{ in.gmtfile | r }}
config <- {{ in.configfile | config: "toml" | r }}
outdir <- {{ out.outdir | r }}
ncores <- {{ envs.ncores | r }}
fgsea <- {{ envs.fgsea | r }}
top <- {{ envs.top | r }}
prerank_method <- {{ envs.prerank_method | r }}

set.seed(8525)
groupby = config$grouping$groupby
if (grepl("^ident", groupby, ignore.case = TRUE)) {
    groupby = "seurat_clusters"
}

sceobj <- readRDS(sceobjfile)

do_one_group <- function(sce, group, outputdir) {
    groupname = if (is.na(as.integer(group))) group else paste0("Cluster", group)
    odir = file.path(outputdir, groupname)
    dir.create(odir, showWarnings = FALSE)
    classes = as.character(sce[[groupby]])
    classes[classes != group] <- "_REST"
    classes[classes == group] <- groupname
    tryCatch({
        if (fgsea) {
            ranks = prerank(
                assay(sce, "exprs"),
                groupname,
                "_REST",
                classes,
                method = prerank_method
            )

            runFGSEA(
                ranks,
                gmtfile,
                top = top,
                outdir = odir
            )
        } else {
            runGSEA(
                assay(sce, "exprs"),
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

do_one_subset <- function(subset) {
    outputdir <- file.path(outdir, subset)
    dir.create(outputdir, showWarnings = FALSE)

    subset_sce <- sceobj[, sceobj$.subset == subset]
    metabolic_sce <- subset_sce[rowData(subset_sce)$metabolic, ]

    groups = metabolic_sce[[groupby]]
    x = mclapply(as.character(unique(groups)), function(group) {
        do_one_group(metabolic_sce, group, outputdir)
    }, mc.cores = ncores)
    if (any(unlist(lapply(x, class)) == "try-error")) {
        stop("mclapply error")
    }
}

subsets <- unique(sceobj$.subset)
for (subset in subsets) {
    do_one_subset(subset)
}
