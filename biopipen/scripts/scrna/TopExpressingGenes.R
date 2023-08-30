library(Seurat)
library(tibble)
library(enrichR)
library(rlang)
library(dplyr)

setEnrichrSite("Enrichr")

srtfile <- {{in.srtobj | r}}
outdir <- {{out.outdir | r}}
mutaters <- {{ envs.mutaters | r }}
ident <- {{ envs.ident | r }}
group.by <- {{ envs["group-by"] | r }}
each <- {{ envs.each | r }}
prefix_each <- {{ envs.prefix_each | r }}
section <- {{ envs.section | r }}
dbs <- {{ envs.dbs | r }}
n <- {{ envs.n | r }}
cases <- {{ envs.cases | r: todot="-" }}

set.seed(8525)

print("- Loading Seurat object ...")
srtobj = readRDS(srtfile)

print("- Mutate meta data if needed ...")
if (!is.null(mutaters) && length(mutaters)) {
    srtobj@meta.data <- srtobj@meta.data %>%
        mutate(!!!lapply(mutaters, parse_expr))
}

print("- Expanding cases ...")
if (is.null(cases) || length(cases) == 0) {
    cases = list(
        Cluster = list(
            ident = ident,
            group.by = group.by,
            each = each,
            prefix_each = prefix_each,
            section = section,
            dbs = dbs,
            n = n
        )
    )
} else {
    cases = lapply(cases, function(cs) {
        if (is.null(cs$ident)) cs$ident = ident
        if (is.null(cs$group.by)) cs$group.by = group.by
        if (is.null(cs$each)) cs$each = each
        if (is.null(cs$prefix_each)) cs$prefix_each = prefix_each
        if (is.null(cs$section)) cs$section = section
        if (is.null(cs$dbs)) cs$dbs = dbs
        if (is.null(cs$n)) cs$n = n
        cs
    })
}

# Expand each and ident
newcases <- list()
for (name in names(cases)) {
    case <- cases[[name]]
    if (is.null(case$each) && !is.null(case$ident)) {
        newcases[[paste0(case$section, ":", name)]] <- case
    } else if (is.null(case$each)) {
        idents <- srtobj@meta.data %>% pull(case$group.by) %>% unique() %>% na.omit()
        for (ident in idents) {
            key = paste0(name, ":", ident)
            newcases[[key]] <- case
            newcases[[key]]$ident <- ident
        }
    } else {
        eachs <- srtobj@meta.data %>% pull(case$each) %>% unique() %>% na.omit()
        for (each in eachs) {
            by <- paste0(".", name, "_", each)
            srtobj@meta.data = srtobj@meta.data %>% mutate(
                !!sym(by) := if_else(
                    !!sym(case$each) == each,
                    !!sym(case$group.by),
                    NA_character_
                )
            )
            case$group.by <- by
            if (is.null(case$ident)) {
                idents <- srtobj@meta.data %>% pull(case$group.by) %>% unique() %>% na.omit()
                for (ident in idents) {
                    key = if (case$prefix_each) {
                        paste0(case$each, "-", each, "-", name, ":", ident)
                    } else {
                        paste0(each, "-", name, ":", ident)
                    }
                    newcases[[key]] <- case
                    newcases[[key]]$ident <- ident
                }
            } else {
                key = if (case$prefix_each) {
                    paste0(case$each, "-", each, ":", name)
                } else {
                    paste0(each, ":", name)
                }
                newcases[[key]] <- case
            }
        }
    }
}
cases <- newcases

do_enrich = function(expr, odir) {
    print("  Saving expressions ...")
    write.table(
        expr %>% as.data.frame() %>% rownames_to_column("Gene"),
        file.path(odir, "expr.txt"),
        sep="\t",
        row.names=TRUE,
        col.names=TRUE,
        quote=FALSE
    )
    write.table(
        expr %>% as.data.frame() %>% rownames_to_column("Gene") %>% head(n),
        file.path(odir, "exprn.txt"),
        sep="\t",
        row.names=TRUE,
        col.names=TRUE,
        quote=FALSE
    )

    print("  Running enrichment ...")
    enriched = enrichr(rownames(head(expr, n)), dbs)
    for (db in dbs) {
        write.table(
            enriched[[db]],
            file.path(odir, paste0("Enrichr-", db, ".txt")),
            sep="\t",
            row.names=FALSE,
            col.names=TRUE,
            quote=FALSE
        )
        png(
            file.path(odir, paste0("Enrichr-", db, ".png")),
            res=100, height=1000, width=1000
        )
        print(plotEnrich(enriched[[db]], showTerms = 20, title=db))
        dev.off()
    }
}

do_case <- function(casename) {
    print(paste("- Running for case:", casename))
    case <- cases[[casename]]

    print("  Calculating average expression ...")
    avgexpr <- AverageExpression(srtobj, group.by = case$group.by)$RNA[, case$ident, drop = FALSE]
    avgexpr <- avgexpr[order(-avgexpr), , drop = FALSE]

    odir = file.path(outdir, casename)
    dir.create(odir, recursive = TRUE, showWarnings = FALSE)

    do_enrich(avgexpr, odir)
}

sapply(names(cases), do_case)
