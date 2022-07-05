source("{{biopipen_dir}}/utils/gsea.R")
library(rlang)
library(Seurat)
library(tidyseurat)

srtfile = {{in.srtobj | r}}
{% if in.casefile %}
cases = {{in.casefile | toml_load | r}}
{% else %}
cases = {{envs | r}}
{% endif %}

outdir = {{out.outdir | r}}
gmtfile = {{envs.gmtfile | r}}
envs = {{envs | r}}

if (is.null(gmtfile) || nchar(gmtfile) == 0) {
    stop("No `envs.gmtfile` provided.")
}

srtobj = readRDS(srtfile)

prepare_exprmat = function(casepms) {
    sobj = srtobj
    if (!is.null(casepms$mutaters)) {
        expr = list()
        for (key in names(casepms$mutaters)) {
            expr[[key]] = parse_expr(casepms$mutaters[[key]])
        }
        sobj = sobj |> mutate(!!!expr)
    }
    if (!is.null(casepms$filter)) {
        sobj = sobj |> filter(eval(parse(text=casepms$filter)))
    }

    samples = rownames(sobj@meta.data[
        sobj@meta.data[[casepms$group.by]] %in% c(casepms$ident.1, casepms$ident.2),
        ,
        drop=FALSE
    ])
    allclasses = sobj@meta.data[samples, casepms$group.by, drop=TRUE]
    exprs = as.data.frame(
        GetAssayData(sobj, slot = "data", assay = "RNA")
    )[, samples, drop=FALSE]
    list(exprs=exprs, allclasses=allclasses)
}

do_case = function(case, casepms) {
    print(paste("- Processing case:", case, "..."))
    odir = file.path(outdir, case)
    dir.create(odir, showWarnings = FALSE)

    exprinfo = prepare_exprmat(casepms)
    ranks = prerank(
        exprinfo$exprs,
        casepms$ident.1,
        casepms$ident.2,
        exprinfo$allclasses,
        envs$method
    )

    write.table(
        ranks,
        file.path(odir, "fgsea.rank"),
        row.names=F,
        col.names=T,
        sep="\t",
        quote=F
    )

    case_envs = envs
    top = case_envs$top
    case_envs$name = NULL
    case_envs$cases = NULL
    case_envs$gmtfile = NULL
    case_envs$nproc = case_envs$ncores
    case_envs$method = NULL
    case_envs$ncores = NULL
    case_envs$top = NULL
    # the rest are the arguments for `fgsea()`
    runFGSEA(ranks, gmtfile, top, odir, case_envs)

}

do_case_with_ident = function(case, casepms, ident) {
    print(paste("- Processing case:", case, "..."))
    print(paste("  Cluster:", ident, "..."))
    ident_name = as.integer(ident)
    if (is.integer(ident_name)) {
        ident_name = paste0("Cluster", ident)
    }
    odir = file.path(outdir, case, ident_name)
    dir.create(odir, showWarnings = FALSE, recursive = TRUE)
    cat("", file = file.path(outdir, case, "percluster"))

    exprinfo = prepare_exprmat(casepms)
    ranks = prerank(
        exprinfo$exprs,
        casepms$ident.1,
        casepms$ident.2,
        exprinfo$allclasses,
        envs$method
    )

    write.table(
        ranks,
        file.path(odir, "fgsea.rank"),
        row.names=F,
        col.names=T,
        sep="\t",
        quote=F
    )

    case_envs = envs
    top = case_envs$top
    case_envs$name = NULL
    case_envs$cases = NULL
    case_envs$gmtfile = NULL
    case_envs$nproc = case_envs$ncores
    case_envs$method = NULL
    case_envs$ncores = NULL
    case_envs$top = NULL
    # the rest are the arguments for `fgsea()`
    runFGSEA(ranks, gmtfile, top, odir, case_envs)

}



.replace_placeholder = function(s, ident) {
    if (is.null(s)) {
        return(s)
    }
    s = sub("{ident}", ident, s, fixed = TRUE)
    s = sub("{cluster}", ident, s, fixed = TRUE)
    s
}

do_case_with_tpl = function(case_with_tpl) {
    casepms = cases$cases[[case_with_tpl]]
    if (isTRUE(casepms$percluster)) {
        # has template in case names
        # currently only cluster is supported
        casepms1 = casepms
        casepms1$percluster = NULL
        for (ident in sort(unique(Idents(srtobj)))) {
            casepms1$ident.1 = .replace_placeholder(casepms$ident.1, ident)
            casepms1$ident.2 = .replace_placeholder(casepms$ident.2, ident)
            casepms1$group.by = .replace_placeholder(casepms$group.by, ident)
            casepms1$filter = .replace_placeholder(casepms$filter, ident)
            if (!is.null(casepms$mutaters)) {
                for (mutname in names(casepms$mutaters)) {
                    casepms1$mutaters[[mutname]] = .replace_placeholder(
                        casepms$mutaters[[mutname]],
                        ident
                    )
                }
            }
            do_case_with_ident(case_with_tpl, casepms1, ident)
        }
    } else {
        do_case(case_with_tpl, casepms)
    }
}

# parallelize?
sapply(names(cases$cases), do_case_with_tpl)
