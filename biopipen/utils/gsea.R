library(ggplot2)
library(dplyr)
library(tibble)

prerank = function(
    exprdata,
    pos,
    neg,
    classes, # must be in the order of colnames(exprdata)
    method = "signal_to_noise"
) {
    library(matrixStats)
    set.seed(8525)
    # See: https://gseapy.readthedocs.io/en/latest/_modules/gseapy/algorithm.html#ranking_metric
    expr_pos_mean = rowMeans(exprdata[, classes == pos, drop=F], na.rm=TRUE)
    expr_neg_mean = rowMeans(exprdata[, classes == neg, drop=F], na.rm=TRUE)
    expr_pos_std = rowSds(as.matrix(exprdata[, classes == pos, drop=F]), na.rm=TRUE)
    expr_neg_std = rowSds(as.matrix(exprdata[, classes == neg, drop=F]), na.rm=TRUE)
    rands = rnorm(length(expr_neg_std)) * 1e-6

    if (method %in% c("s2n", "signal_to_noise")) {
        out = (expr_pos_mean - expr_neg_mean) / (expr_pos_std + expr_neg_std + rands)
    } else if (method %in% c("abs_s2n", "abs_signal_to_noise")) {
        out = abs((expr_pos_mean - expr_neg_mean) / (expr_pos_std + expr_neg_std + rands))
    } else if (method == "t_test") {
        # ser = (df_mean[pos] - df_mean[neg])/ np.sqrt(df_std[pos]**2/len(df_std)+df_std[neg]**2/len(df_std) )
        out = (expr_pos_mean - expr_neg_mean) / sqrt(
            expr_pos_std ^ 2 / length(expr_pos_std) +
            expr_neg_std ^ 2 / length(expr_neg_std)
        )
    } else if (method == "ratio_of_classes") {
        out = expr_pos_mean / expr_neg_mean
    } else if (method == "diff_of_classes") {
        out = expr_pos_mean - expr_neg_mean
    } else if (method == "log2_ratio_of_classes") {
        out = log2(expr_pos_mean) - log2(expr_neg_mean)
    } else {
        stop(paste("Unknown method:", method))
    }
    # todo: log2fc * -log10(p)
    # see https://github.com/crazyhottommy/RNA-seq-analysis/blob/master/GSEA_explained.md#2-using-a-pre-ranked-gene-list
    out = as.data.frame(out) %>% rownames_to_column("Gene") %>% arrange(.[[2]])
    colnames(out)[2] = paste(pos, "vs", neg, sep="_")
    return(out)
}

runEnrichr = function(
    genes,
    dbs,
    outdir,
    showTerms = 20,
    numChar =40,
    orderBy = "P.value"
) {
    library(enrichR)
    setEnrichrSite("Enrichr") # Human genes

    enriched = enrichr(genes, dbs)

    for (db in dbs) {
        enr = enriched[[db]] %>% select(-c(Old.P.value, Old.Adjusted.P.value))
        outtable = file.path(outdir, paste0("Enrichr_", db, ".txt"))
        outfig = file.path(outdir, paste0("Enrichr_", db, ".png"))
        write.table(enr, outtable, row.names=T, col.names=F, sep="\t", quote=F)

        png(outfig, res=100, height=1000, width=1400)
        print(
            plotEnrich(
                enriched[[db]],
                showTerms=showTerms,
                numChar=numChar,
                orderBy=orderBy
            )
        )
        dev.off()
    }

}

runFGSEA = function(
    ranks,
    gmtfile,
    top,
    outdir,
    envs = list(),
    plot = TRUE  # only generate fgsea.txt?
) {
    library(data.table)
    library(fgsea)
    set.seed(8525)

    if (is.data.frame(ranks)) {
        ranks = setNames(ranks[[2]], ranks[[1]])
    } else if (is.list(ranks)) {
        ranks = unlist(ranks)
    }

    envs$pathways = gmtPathways(gmtfile)
    envs$stats = ranks
    gsea_res = do.call(fgsea::fgsea, envs)
    gsea_res = gsea_res[order(pval), ]

    write.table(
        gsea_res %>%
            mutate(leadingEdge = sapply(leadingEdge, function(x) paste(x, collapse=","))),
        file = file.path(outdir, "fgsea.txt"),
        row.names = FALSE,
        col.names = TRUE,
        sep = "\t",
        quote = FALSE
    )

    if (!plot) {return (NULL)}

    if (top > 1) {
        topPathways = head(gsea_res, n=top)[, "pathway"]
    } else {
        topPathways = gsea_res[padj < top][, "pathway"]
    }
    topPathways = unlist(topPathways)

    tablefig = file.path(outdir, "gsea_table.png")
    png(tablefig, res=100, width=1000, height=200 + 40 * length(topPathways))
    print(plotGseaTable(
        envs$pathways[topPathways],
        ranks,
        gsea_res,
        gseaParam = if (!is.null(envs$gseaParam)) envs$gseaParam else 1
    ))
    dev.off()

    for (pathway in topPathways) {
        enrfig = file.path(outdir, paste0("fgsea_", gsub("/", "-", pathway, fixed=T), ".png"))
        png(enrfig, res=100, width=1000, height=800)
        print(plotEnrichment(
            envs$pathways[[pathway]],
            ranks,
            gseaParam = if (!is.null(envs$gseaParam)) envs$gseaParam else 1
        ) + labs(title = pathway))
        dev.off()
    }
}

runGSEA = function(
    indata, # expression data
    classes, # sample classes
    gmtfile, # the GMT file
    outdir,
    envs = list() # other arguments for GSEA()
) {
    library(GSEA)
    # reproducibility
    if (is.null(envs$random.seed)) {
        envs$random.seed <- 8525
    }

    # prepare gct file
    gctfile = file.path(outdir, "gsea.gct")
    con = file(gctfile, open='w')
    write("#1.2", con)
    write(paste(dim(indata), collapse = "\t"), con)
    close(con)
    indata = indata %>%
        as.data.frame() %>%
        mutate(Description = "na") %>%
        rownames_to_column("NAME") %>%
        select(NAME, Description, everything())
    write.table(
        indata,
        gctfile,
        row.names = F,
        col.names = T,
        sep="\t",
        quote=F,
        append = T
    )

    # prepare cls file
    clsfile = file.path(outdir, "gsea.cls")
    uniclasses = unique(classes)
    con = file(clsfile, open='w')
    write(paste(length(classes), length(uniclasses), '1'), con)
    write(paste('#', paste(uniclasses, collapse=" ")), con)
    write(paste(classes, collapse=" "), con)
    close(con)

    envs$input.ds = gctfile
    envs$input.cls = clsfile
    envs$gs.db = gmtfile
    envs$output.directory = outdir

    do.call(GSEA, envs)
}
