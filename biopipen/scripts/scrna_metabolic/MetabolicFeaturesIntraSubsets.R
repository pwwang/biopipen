source("{{biopipen_dir}}/utils/gsea.R")

library(parallel)
library(scater)

sceobjfile <- {{ in.sceobj | r }}
gmtfile <- {{ in.gmtfile | r }}
config <- {{ in.configfile | read | toml_loads | r }}
outdir <- {{ out.outdir | r }}
ncores <- {{ envs.ncores | r }}
fgsea <- {{ envs.fgsea | r }}
top <- {{ envs.top | r }}
prerank_method <- {{ envs.prerank_method | r }}

set.seed(8525)

sceobj <- readRDS(sceobjfile)

do_one_design <- function(sce, dsname, case, control, groupdir) {
    print(paste("  Design:", dsname, "(", case, ",", control, ")"))
    sce_case = sce[, sce$.subset == case]
    sce_control = sce[, sce$.subset == control]
    exprs_case = assay(sce_case, "exprs")
    exprs_control = assay(sce_control, "exprs")
    if (ncol(exprs_case) < 3 || ncol(exprs_control) < 3) {
        print("          Skip (not enough cells)")
        return (NULL)
    }
    odir = file.path(groupdir, dsname)
    dir.create(odir, showWarnings = FALSE)
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
            outdir = odir
        )
    } else {
        runGSEA(
            cbind(exprs_case, exprs_control),
            c(rep(case, ncol(exprs_case)), rep(control, ncol(exprs_control))),
            gmtfile,
            odir
        )
    }
}

do_one_group <- function(group) {
    print(paste("- Group:", group, "..."))
    sce = sceobj[rowData(sceobj)$metabolic, as.character(sceobj[[config$grouping_name]]) == group]
    groupname = if (is.na(as.integer(group))) group else paste0("Cluster", group)
    groupdir = file.path(outdir, groupname)
    dir.create(groupdir, showWarnings = FALSE)

    for (dsname in names(config$design)) {
        do_one_design(
            sce,
            dsname,
            config$design[[dsname]][1],
            config$design[[dsname]][2],
            groupdir
        )
    }
}

mclapply(
    as.character(unique(sceobj[[config$grouping_name]])),
    do_one_group,
    mc.cores = ncores
)
