{{ "plot.R", "__init__.R" | rimport }}

infile <- {{ i.infile | quote }}
outfile <- {{ o.outfile | quote }}
fccut <- {{ args.fccut | R }}
pcut <- {{ args.pcut | R }}
usep <- {{ args.usepval | R }}
hilights <- {{ args.hilights | R }}
devpars <- {{ args.devpars | R }}
ggs <- {{ args.ggs | R }}
fdrpos <- {{ args.fdrpos | R }}

indata <- read.table.inopts(infile,
                            list(rnames = TRUE, cnames = TRUE, dup = "ignore"))
plotdata <- data.frame(.placeholder = 1:nrow(indata))
rownames(plotdata) <- rownames(indata)
for (cname in colnames(indata)) {
    cname_to_check <- tolower(gsub("[^[:alnum:] ]", "", cname))
    if (cname_to_check %in% c("logfc", "log2fc",
                              "logfoldchange", "log2foldchange")) {
        plotdata$logFC = indata[, cname]
    } else if (usep && (cname_to_check %in% c("p", "pvalue", "pval"))) {
        plotdata$PVal = indata[, cname]
    } else if (!usep && (cname_to_check %in% c("q", "qvalue", "qval", "fdr",
                                               "adjpval", "adjpvalue",
                                               "padj"))) {
        plotdata$FDR = indata[, cname]
    } else {
        #plotdata[[cname]] = indata[, cname]
    }
}

plot.volplot(plotdata[, -1, drop=FALSE], outfile,
             params = list(logfccut=fccut, pcut=pcut, hilight=hilights,
                           fdrpos = fdrpos),
             ggs = ggs, devpars = devpars)
