source("{{biopipen_dir}}/utils/misc.R")

library(rlang)
library(dplyr)
library(tidyr)
library(tibble)
library(Seurat)
library(enrichR)
library(ggplot2)
library(future)
library(tidyseurat)

setEnrichrSite("Enrichr")

srtobjfile = {{in.srtobj | quote}}
outdir = {{out.outdir | quote}}
{% if in.casefile %}
cases = {{in.casefile | toml_load | r: todot="-"}}
{% else %}
cases = {{envs | r: todot="-"}}
{% endif %}
dbs = {{envs.dbs | r}}
ncores = {{envs.ncores | r}}
sigmarkers = {{envs.sigmarkers | r}}

set.seed(8525)
options(future.globals.maxSize = 80000 * 1024^2)
plan(strategy = "multicore", workers = ncores)

if (length(cases) == 0) {
    stop("No `envs.cases` or `in.casefile` provided.")
}

seurat_obj = readRDS(srtobjfile)

do_enrich = function(case, markers) {
    print(paste("  Running enrichment for case:", case))
    casedir = file.path(outdir, case)
    dir.create(casedir, showWarnings = FALSE)
    if (nrow(markers) == 0) {
        print(paste("  No markers found for case:", case))
        cat("No markers found.", file=file.path(casedir, "error.txt"))
        return()
    }
    markers_sig = markers %>% filter(!!parse_expr(sigmarkers))
    if (nrow(markers_sig) == 0) {
        print(paste("  No significant markers found for case:", case))
        cat("No significant markers.", file=file.path(casedir, "error.txt"))
        return()
    }
    write.table(
        markers_sig,
        file.path(casedir, "markers.txt"),
        sep="\t",
        row.names=FALSE,
        col.names=TRUE,
        quote=FALSE
    )
    if (nrow(markers_sig) < 5) {
        for (db in dbs) {
            write.table(
                data.frame(Warning = "Not enough significant markers."),
                file.path(casedir, paste0("Enrichr-", db, ".txt")),
                sep="\t",
                row.names=FALSE,
                col.names=TRUE,
                quote=FALSE
            )
            png(
                file.path(casedir, paste0("Enrichr-", db, ".png")),
                res=100, height=200, width=1000
            )
            print(
                ggplot() +
                annotate("text", x=1, y=1, label="Not enough significant markers.") +
                theme_classic()
            )
            dev.off()
        }
    } else {
        enriched = enrichr(markers_sig$gene, dbs)
        for (db in dbs) {
            write.table(
                enriched[[db]],
                file.path(casedir, paste0("Enrichr-", db, ".txt")),
                sep="\t",
                row.names=FALSE,
                col.names=TRUE,
                quote=FALSE
            )
            png(
                file.path(casedir, paste0("Enrichr-", db, ".png")),
                res=100, height=1000, width=1000
            )
            print(plotEnrich(enriched[[db]], showTerms = 20, title=db))
            dev.off()
        }
    }
}

mutate_meta = function(obj, mutaters) {
    meta = obj@meta.data
    if (!is.null(mutaters)) {
        expr = list()
        for (key in names(mutaters)) {
            expr[[key]] = parse_expr(mutaters[[key]])
        }
        obj@meta.data = meta %>% mutate(!!!expr)
    }
    return(obj)
}

do_case = function(case) {
    cat(paste("- Dealing with case:", case, "...\n"))
    casepms = cases$cases[[case]]
    obj = seurat_obj
    if (!is.null(casepms$filter)) {
        obj = obj %>% filter(eval(parse(text=casepms$filter)))
    }
    obj = mutate_meta(obj, casepms$mutaters)
    casepms$mutaters = NULL
    if (!is.null(casepms$filter2)) {
        obj = obj %>% filter(eval(parse(text=casepms$filter2)))
    }
    if (!is.null(casepms$each)) {
        eachs = unique(obj@meta.data[[casepms$each]])
        for (each in eachs) {
            print(paste("  Dealing with unit:", each, "..."))
            eachobj = obj %>% filter(!!parse_expr(casepms$each) == each)
            casepms$object = eachobj
            markers = do_call(FindMarkers, casepms) %>% rownames_to_column("gene")
            do_enrich(paste0(case, " (", each, ")"), markers)
        }
    } else {
        if (is.null(casepms$ident.1) && is.null(casepms$ident.2)) {
            Idents(obj) = casepms$group.by
            casepms$group.by = NULL
            casepms$object = obj
            allmarkers = do_call(FindAllMarkers, casepms)
            # Is it always cluster?
            for (group in sort(unique(allmarkers$cluster))) {
                do_enrich(paste(case, group, sep="_"), allmarkers %>% filter(cluster == group))
            }
        } else {
            casepms$object = obj
            markers = do_call(FindMarkers, casepms) %>% rownames_to_column("gene")
            do_enrich(case, markers)
        }
    }
}

sapply(names(cases$cases), do_case)
