{{ biopipen_dir | joinpaths: "utils", "misc.R" | source_r }}

library(rlang)
library(stats)
library(ggplot2)
library(ggprism)
library(qqplotr)

theme_set(theme_prism())

infile <- {{in.infile | r}}
theorfile <- {{in.theorfile | r}}
outfile <- {{out.outfile | r}}
val_col <- {{envs.val_col | r}}
theor_col <- {{envs.theor_col | r}}
theor_trans <- {{envs.theor_trans | r}}
theor_funs <- {{envs.theor_funs | r}}
devpars <- {{envs.devpars | r}}
title <- {{envs.title | r}}
xlabel <- {{envs.xlabel | r}}
ylabel <- {{envs.ylabel | r}}
kind <- {{envs.kind | r}}
trans <- {{envs.trans | r}}
args <- {{envs.args | r}}
band_args <- {{envs.band | r}}
line_args <- {{envs.line | r}}
point_args <- {{envs.point | r}}
ggs <- {{envs.ggs | r}}

.eval_fun <- function(fun) {
    if (is.character(fun)) {
        fun <- trimws(fun)
        if (grepl("^-\\s*[a-zA-Z\\.][0-9a-zA-Z\\._]*$", fun)) {
            fun <- trimws(substring(fun, 2))
            fun <- eval(parse(text = fun))
            return(function(x) -fun(x))
        } else {
            return(eval(parse(text = fun)))
        }
    } else {
        return(fun)
    }
}

indata <- read.table(infile, header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names = FALSE)
if (is.numeric(val_col)) {
    val_col <- colnames(indata)[val_col]
}
if (!is.null(trans)) {
    trans <- .eval_fun(trans)
    indata[[val_col]] <- trans(indata[[val_col]])
}

if (!is.null(theor_col)) {
    if (is.numeric(theor_col)) {
        theor_col <- colnames(theor)[theor_col]
    }

    if (!is.null(theorfile)) {
        theor <- read.table(theorfile, header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names = FALSE)
        theor_vals <- theor[[theor_col]]
    } else {
        theor_vals <- indata[[theor_col]]
    }

    if (!is.null(theor_trans)) {
        theor_trans <- .eval_fun(theor_trans)
        theor_vals <- theor_trans(theor_vals)
    }
    theor_vals <- sort(na.omit(theor_vals))
}

band_fun <- ifelse(kind == "pp", stat_pp_band, stat_qq_band)
line_fun <- ifelse(kind == "pp", stat_pp_line, stat_qq_line)
point_fun <- ifelse(kind == "pp", stat_pp_point, stat_qq_point)

for (fun in names(theor_funs)) {
    assign(fun, .eval_fun(theor_funs[[fun]]))
}

if (!is.null(band_args) || isFALSE(band_args)) {
    if (isTRUE(band_args$disabled)) {
        band_args <- NULL
    } else {
        band_args$disabled <- NULL
        band_args <- list_update(band_args, args)
        if (band_args$distribution == "custom") {
            band_args$dparams <- band_args$dparams %||% list()
            band_args$dparams$values <- theor_vals
        }
    }
}
if (!is.null(line_args) || isFALSE(line_args)) {
    if (isTRUE(line_args$disabled)) {
        line_args <- NULL
    } else {
        line_args$disabled <- NULL
        line_args <- list_update(line_args, args)
        if (line_args$distribution == "custom") {
            line_args$dparams <- line_args$dparams %||% list()
            line_args$dparams$values <- theor_vals
        }
    }
}
if (!is.null(point_args) || isFALSE(point_args)) {
    if (isTRUE(point_args$disabled)) {
        point_args <- NULL
    } else {
        point_args$disabled <- NULL
        point_args <- list_update(point_args, args)
        if (point_args$distribution == "custom") {
            point_args$dparams <- point_args$dparams %||% list()
            point_args$dparams$values <- theor_vals
        }
    }
}

title <- title %||% waiver()
xlabel <- xlabel %||% waiver()
ylabel <- ylabel %||% waiver()

indata <- indata[complete.cases(indata), , drop = FALSE]
indata <- indata[order(indata[[val_col]]), , drop = FALSE]

p <- ggplot(data = indata, mapping = aes(sample = !!sym(val_col))) +
    labs(title = title, x = xlabel, y = ylabel)

if (!is.null(band_args)) {
    p <- p + do_call(band_fun, band_args)
}
if (!is.null(line_args)) {
    p <- p + do_call(line_fun, line_args)
}
if (!is.null(point_args)) {
    p <- p + do_call(point_fun, point_args)
}

if (!is.null(ggs)) {
    for (gg in ggs) {
        p <- p + eval(parse(text = gg))
    }
}

png(outfile, width=devpars$width, height=devpars$height, res=devpars$res)
print(p)
dev.off()
