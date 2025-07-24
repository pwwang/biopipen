library(biopipen.utils)

infiles <- {{in.infiles | r}}
outfile <- {{out.outfile | r}}
sep <- {{envs.sep | r}}
header <- {{envs.header | r}}
filenames <- {{envs.filenames | r}}
filenames_col <- {{envs.filenames_col | r}}

filenames_fns <- c()
if (!is.null(filenames)) {
    if (!is.character(filenames)) {
        stop("`envs.filenames` must be a character vector or a string of R `function`")
    }
    filenames <- trimws(filenames)
    if (length(filenames) == 1) {
        if (startsWith(filenames, "function")) {
            filenames_fns <- eval(parse(text = filenames))
        } else {
            filenames_fns <- trimws(strsplit(filenames, ",")[[1]])
        }
    } else {
        filenames_fns <- filenames
    }
    if (!is.function(filenames_fns) && length(filenames_fns) < length(infiles)) {
        stop("Length of `envs.filenames` must be equal to or greater than the number of input files")
    }
}

indata <- list()
for (i in seq_along(infiles)) {
    indata[[i]] <- read.delim(infiles[[i]], sep = sep, header = header)
    if (is.null(filenames)) {
        next
    }
    if (is.function(filenames_fns)) {
        indata[[i]][filenames_col] <- filenames_fns(infiles[i])
    } else {
        indata[[i]][filenames_col] <- filenames_fns[i]
    }
}

outdata = do_call(rbind, indata)

write.table(
    outdata,
    file = outfile,
    sep = sep,
    row.names = FALSE,
    col.names = header,
    quote = FALSE
)
