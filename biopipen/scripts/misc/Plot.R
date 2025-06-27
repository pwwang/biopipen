library(gglogger)
library(plotthis)
library(rlang)
library(biopipen.utils)

datafile <- {{in.datafile | r}}
plotfile <- {{out.plotfile | r}}
plotprefix <- {{out.plotfile | prefix | r}}
read_opts <- {{envs.read_opts | r: todot="-"}}
envs <- {{envs | r}}

fn <- envs$fn
envs$fn <- NULL
devpars <- envs$devpars
envs$devpars <- NULL
more_formats <- envs$more_formats
envs$more_formats <- NULL
save_code <- envs$save_code
envs$save_code <- NULL
envs$read_opts <- NULL

if (endsWith(datafile, ".qs") || endsWith(datafile, ".qs2") ||
    endsWith(datafile, ".rds") || endsWith(datafile, ".RDS")) {
    envs$data <- read_obj(datafile)
} else {
    read_opts <- read_opts %||% list()
    read_opts$file <- datafile
    envs$data <- do.call(read.table, read_opts)
}

if (fn == "ManhattanPlot" && !is.null(envs$chromosomes)) {
    norm_chroms <- function(chrs) {
        chrs <- as.character(chrs)
        if (length(chrs) == 1 && grepl(",", chrs)) {
            chrs <- trimws(unlist(strsplit(chrs, ",")))
        }
        if (length(chrs) > 1) {
            return(unique(unlist(sapply(chrs, function(chr) norm_chroms(chr)))))
        }
        if (!grepl("-", chrs)) { return(chrs) }

        # expand chr1-22 -> chr1, chr2, ..., chr22
        # chr1-22 -> 'chr1', '22'
        chrs <- unlist(strsplit(chrs, "-"))
        if (length(chrs) != 2) {
            stop(paste0("Invalid chroms: ", chrs))
        }
        # detect prefix
        prefix1 <- gsub("[0-9]", "", chrs[1])
        prefix2 <- gsub("[0-9]", "", chrs[2])
        if (nchar(prefix2) > 0 && prefix1 != prefix2) {
            stop(paste0("Invalid chroms: ", chrs, " (prefix mismatch)"))
        }
        chr_a <- as.integer(substring(chrs[1], nchar(prefix1) + 1))
        chr_b <- as.integer(substring(chrs[2], nchar(prefix2) + 1))
        chr_min <- min(chr_a, chr_b)
        chr_max <- max(chr_a, chr_b)
        return(paste0(prefix1, chr_min:chr_max))
    }

    envs$chromosomes <- norm_chroms(envs$chromosomes)
}

plotfn <- utils::getFromNamespace(fn, "plotthis")
if (save_code) {
    plotfn <- gglogger::register(plotfn, name = fn)
}

p <- do_call(plotfn, envs)
save_plot(p, plotprefix, devpars, formats = unique(c("png", more_formats)))

if (save_code) {
    save_plotcode(
        p,
        setup = c('library(plotthis)', '', 'load("data.RData")', 'list2env(envs, envir = .GlobalEnv)'),
        prefix = plotprefix,
        "envs",
        auto_data_setup = FALSE
    )
}
