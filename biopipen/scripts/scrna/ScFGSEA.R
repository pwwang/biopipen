source("{{biopipen_dir}}/utils/gsea.R")
library(rlang)
library(Seurat)

srtfile = {{in.srtobj | r}}
{% if in.casefile %}
cases = {{in.casefile | toml_load | r}}
{% else %}
cases = {{envs.cases | r}}
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

do_case = function(case) {
    odir = file.path(outdir, case)
    dir.create(odir, showWarnings = FALSE)
    casepms = cases$cases[[case]]
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
    case_envs$cases = NULL
    case_envs$gmtfile = NULL
    case_envs$nproc = case_envs$ncores
    case_envs$method = NULL
    case_envs$ncores = NULL
    case_envs$top = NULL
    # the rest are the arguments for `fgsea()`
    runFGSEA(ranks, gmtfile, top, odir, case_envs)

}

# parallelize?
sapply(names(cases$cases), do_case)
