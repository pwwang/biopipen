library(Seurat)
library(tibble)
library(enrichR)

setEnrichrSite("Enrichr")

srtfile = {{in.srtobj | r}}
outdir = {{out.outdir | r}}
cluster_col = {{envs.cluster_col | r}}
assay = {{envs.assay | r}}
dbs = {{envs.dbs | r}}
n = {{envs.n | r}}

print("- Loading Seurat object ...")
srtobj = readRDS(srtfile)

print("- Calculating average expression ...")
Idents(srtobj) = cluster_col
avgexpr = AverageExpression(srtobj, assays = assay)[[assay]]
idents = unique(Idents(srtobj))


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


run_ident = function(ident) {
    print(paste("- Running for ident:", ident))
    identname = if (is.na(is.numeric(ident))) ident else paste0("Cluster", ident)
    odir = file.path(outdir, identname)
    dir.create(odir, showWarnings = FALSE)

    expr = avgexpr[, ident, drop = FALSE]
    expr = expr[order(-expr), , drop = FALSE]

    do_enrich(expr, odir)
}

sapply(idents, run_ident)
