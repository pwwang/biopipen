{{ "__init__.R" | rimport }}

infiles <- {{ in.infiles | R }}
outfile <- {{ out.outfile | R }}
na <- {{ args.na | R }}
fn2rname <- {{ args.fn2rname | str | .startswith: "function" ? !R }}
inopts <- {{ args.inopts | R }}
fill <- as.logical({{ args.fill | R }})

mats <- list()
for (i in 1:length(infiles)) {
    mats[[i]] <- read.table.inopts(infiles[i], inopts)

    if (is.function(fn2rname)) {
        rnames <- rownames(mats[[i]])
        fn <- tools::file_path_sans_ext(basename(infiles[i]))
        if (length(formals(fn2rname)) == 1) {
            fnc <- fn2rname(fn)
        } else {
            fnc <- fn2rname(fn, rnames)
        }
        rownames(mats[[i]]) <- fnc
    }
}
mat <- ifelse(fill, do.call(rbind.fill, mats), do.call(rbind, mats))
mat[is.na(mat)] <- na

write.table(mat, outfile, sep = "\t", quote = FALSE,
            row.names = is.function(fn2rname), col.names = inopts$cnames)
