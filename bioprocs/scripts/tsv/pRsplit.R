{{ "__init__.R" | rimport }}

infile <- {{ i.infile | R }}
outdir <- {{ o.outdir | R }}
inopts <- {{ args.inopts | R }}
size <- {{ args.size | int | R }}
ext <- {{ i.infile | ext | R }}
stem <- {{ i.infile | stem | R }}

mat <- read.table.inopts(infile, inopts)
index <- 0

for (chunk in split(1:nrow(mat),
                    ceiling(seq_along(1:nrow(mat)) / size))) {
    index <- index + 1
    # // rnames <- allrnames[chunk]
    m <- mat[chunk, , drop = F]
    # // fn <- paste(gsub("[[:punct:]]", "_", rnames), collapse = "-")
    fn <- paste0(stem, "_ROWS", index, ext)
    outfile <- file.path(outdir, fn)
    outparams <- list(
        x = m,
        file = outfile,
        sep = list.get(inopts, "delimit", "\t"),
        quote = F,
        col.names = is.true(inopts$cnames),
        row.names = is.true(inopts$rnames)
    )
    do.call(write.table, outparams)
}
