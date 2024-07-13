{{ biopipen_dir | joinpaths: "utils", "misc.R" | source_r }}

library(ggpmisc)
library(rlang)
library(ggplot2)
library(ggprism)

theme_set(theme_prism())

infile <- {{in.infile | r}}
outfile <- {{out.outfile | r}}
x_col <- {{envs.x_col | r}}
y_col <- {{envs.y_col | r}}
devpars <- {{envs.devpars | r}}
args <- {{envs.args | r}}
ggs <- {{envs.ggs | r}}
formula <- {{envs.formula | r}}
mapping <- {{envs.mapping | r}}
stats <- {{envs.stats | r}}

.ensure_r <- function(ex, recursive=TRUE) {
    if (is.character(ex)) {
        ex <- trimws(ex)
        if (grepl("^-\\s*[a-zA-Z\\.][0-9a-zA-Z\\._]*$", ex)) {
            ex <- trimws(substring(ex, 2))
            ex <- eval(parse(text = ex))
            return(function(x) -ex(x))
        } else {
            return(eval(parse(text = ex)))
        }
    } else if (is.list(ex) && recursive) {
        return(lapply(ex, .ensure_r, recursive=TRUE))
    } else {
        return(ex)
    }
}

.merge_aes <- function(aes1, aes2) {
    if (is.null(aes1)) {
        return(aes2)
    }
    if (is.null(aes2)) {
        return(aes1)
    }
    merged <- c(aes1, aes2)  # list
    out <- list()
    for (key in names(merged)) {
        if (is.null(out[[key]])) {
            out[[key]] <- merged[[key]]
        } else {
            log_warn(paste("Overwriting mapping key:", key))
        }
    }
    return(do.call(aes, out))
}

if (is.null(formula)) {
    stop("Formula must be provided")
}
if (!is.null(mapping)) {
    if (startsWith(mapping, "(") && endsWith(mapping, ")")) {
        mapping <- paste0("aes", mapping)
    } else if (!startsWith(mapping, "aes(")) {
        mapping <- paste0("aes(", mapping, ")")
    }
    mapping <- .ensure_r(mapping)
}
formula <- as.formula(formula)

indata <- read.table(infile, header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names = FALSE)
if (is.numeric(x_col)) {
    x_col <- colnames(indata)[x_col]
}
if (is.numeric(y_col)) {
    y_col <- colnames(indata)[y_col]
}

args <- lapply(args, .ensure_r)
args$mapping <- .merge_aes(args$mapping, mapping)

if (!is.null(stats)) {
    stats <- lapply(stats, .ensure_r)
}

p <- ggplot(indata, aes(x = !!sym(x_col), y = !!sym(y_col))) +
    do.call(geom_point, args)

for (stat in names(stats)) {
    if (startsWith(stat, "stat_")) {
        stat <- substring(stat, 6)
    }
    if (grepl("#", stat)) {
        st <- paste0("stat_", strsplit(stat, "#")[[1]][1])
    } else {
        st <- paste0("stat_", stat)
    }
    stats[[stat]]$formula <- stats[[stat]]$formula %||% formula
    stats[[stat]]$mapping <- .merge_aes(stats[[stat]]$mapping, mapping)
    p <- p + do.call(st, stats[[stat]])
}

if (!is.null(ggs)) {
    for (gg in ggs) {
        p <- p + eval(parse(text = gg))
    }
}

p <- p + scale_color_biopipen()

png(outfile, width=devpars$width, height=devpars$height, res=devpars$res)
print(p)
dev.off()
