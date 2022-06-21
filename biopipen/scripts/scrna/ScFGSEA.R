source("{{biopipen_dir}}/utils/gsea.R")
library(rlang)
library(Seurat)

srtfile = {{in.srtobj | r}}
{% if in.casefile %}
cases = {{in.casefile | toml_load | r}}
{% else %}
cases = {{envs | r}}
{% endif %}

outdir = {{out.outdir | r}}
gmtfile = {{envs.gmtfile | r}}
envs = {{envs | r}}

srtobj = readRDS(srtfile)

prepare_exprmat = function(casepms) {
    if (!is.null(casepms$mutaters)) {
        expr = list()
        for (key in names(casepms$mutaters)) {
            expr[[key]] = parse_expr(casepms$mutaters[[key]])
        }
        metadata = srtobj@meta.data |> mutate(!!!expr)
    } else {
        metadata = srtobj@meta.data
    }
    samples = rownames(metadata[
        metadata[[casepms$group.by]] %in% c(casepms$ident.1, casepms$ident.2),
        ,
        drop=FALSE
    ])
    allclasses = metadata[samples, casepms$group.by, drop=TRUE]
    exprs = as.data.frame(
        GetAssayData(srtobj, slot = "data", assay = "RNA")
    )[, samples, drop=FALSE]
    list(exprs=exprs, allclasses=allclasses)
}

do_case = function(case, casepms) {
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

.replace_placeholder = function(s, ident) {
    s = sub("{ident}", ident, s, fixed = TRUE)
    s = sub("{cluster}", ident, s, fixed = TRUE)
    s
}

do_case_with_tpl = function(case_with_tpl) {
    if (grepl("{", case_with_tpl, fixed = TRUE)) {
        # has template in case names
        # currently only cluster is supported
        casepms = cases$cases[[case_with_tpl]]
        for (ident in unique(Idents(srtobj))) {
            case = .replace_placeholder(case_with_tpl, ident)
            casepms$ident.1 = .replace_placeholder(casepms$ident.1, ident)
            casepms$ident.2 = .replace_placeholder(casepms$ident.2, ident)
            casepms$group.by = .replace_placeholder(casepms$group.by, ident)
            if (!is.null(casepms$mutaters)) {
                for (mutname in names(casepms$mutaters)) {
                    casepms$mutaters[[mutname]] = .replace_placeholder(
                        casepms$mutaters[[mutname]],
                        ident
                    )
                }
            }
            do_case(case, casepms)
        }
    } else {
        do_case(case_with_tpl, casepms)
    }
}

# parallelize?
sapply(names(cases$cases), do_case_with_tpl)
