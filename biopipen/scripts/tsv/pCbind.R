{{ "__init__.r" | rimport }}

infiles <- {{ i.infiles | R }}
outfile <- {{ o.outfile | R }}
na <- {{ args.na | R }}
fn2cname <- {{ args.fn2cname | str | ?.startswith: "function" | !R }}
inopts <- {{ args.inopts | R }}
fill <- as.logical({{ args.fill | R }})

mats <- list()
for (i in 1:length(infiles)) {
    mats[[i]] <- read.table.inopts(infiles[i], inopts)

    if (is.function(fn2cname)) {
        cnames <- colnames(mats[[i]])
        fn <- tools::file_path_sans_ext(basename(infiles[i]))
        if (length(formals(fn2cname)) == 1) {
            fnc <- fn2cname(fn)
        } else {
            fnc <- fn2cname(fn, cnames)
        }
        colnames(mats[[i]]) <- fnc
    }
}
mat <- ifelse(fill, do.call(cbind.fill, mats), do.call(cbind, mats))
mat[is.na(mat)] <- na

write.table(mat, outfile, sep = "\t", quote = FALSE,
            col.names = is.function(fn2cname), row.names = inopts$rnames)
