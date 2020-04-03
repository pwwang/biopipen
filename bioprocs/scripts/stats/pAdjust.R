
{{ "__init__.R" | rimport }}
infiles <- {{ i.infiles | R }}
outfile <- {{ o.outfile | quote }}
method <- {{ args.method | quote }}
pcol <- {{ args.pcol | R }}
inopts <- {{ args.inopts | R }}

pvals <- c()
for (infile in infiles) {
    indata <- read.table.inopts(infile, inopts)
    pvals <- c(pvals, unlist(indata[, pcol, drop = TRUE]))
}

padjs <- p.adjust(pvals, method)

prev <- 1
for (infile in infiles) {
    indata <- read.table.inopts(infile, inopts)
    indata$Padj <- padjs[prev:(prev + nrow(indata) - 1)]
    write.table(
        indata, outfile,
        append = TRUE, quote = FALSE, sep = "\t",
        col.names = prev == 1 && is.true(inopts$cnames),
        row.names = is.true(inopts$rnames)
    )
    prev <- prev + nrow(indata)
}
