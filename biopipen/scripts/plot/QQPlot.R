source("{{biopipen_dir}}/utils/misc.R")

library(ggplot2)
library(ggprism)
library(qqplotr)

theme_set(theme_prism())

infile <- {{in.infile | r}}
outfile <- {{out.outfile | r}}
val_col <- {{envs.val_col | r}}
devpars <- {{envs.devpars | r}}
title <- {{envs.title | r}}
xlabel <- {{envs.xlabel | r}}
ylabel <- {{envs.ylabel | r}}
kind <- {{envs.kind | r}}
trans <- {{envs.trans | r}}
band_args <- {{envs.band | r}}
line_args <- {{envs.line | r}}
point_args <- {{envs.point | r}}
ggs <- {{envs.ggs | r}}

indata <- read.table(infile, header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names = FALSE)
if (is.numeric(val_col)) { val_col <- colnames(indata)[val_col] }

band_fun <- ifelse(kind == "pp", stat_pp_band, stat_qq_band)
line_fun <- ifelse(kind == "pp", stat_pp_line, stat_qq_line)
point_fun <- ifelse(kind == "pp", stat_pp_point, stat_qq_point)

title <- title %||% waiver()
xlabel <- xlabel %||% waiver()
ylabel <- ylabel %||% waiver()

if (!is.null(trans)) {
    trans <- trimws(trans)
    if (trans == "-log10") {
        trans <- function(x) -log10(x)
    } else {
        trans <- eval(parse(text = trans))
    }

    indata$.trans_val <- trans(indata[[val_col]])
    val_col <- ".trans_val"
}

indata <- indata[!is.na(indata[[val_col]]), , drop=FALSE]

p <- ggplot(data = indata, mapping = aes(sample = !!sym(val_col))) +
    do_call(band_fun, band_args) +
    do_call(line_fun, line_args) +
    do_call(point_fun, point_args) +
    labs(title = title, x = xlabel, y = ylabel)

if (!is.null(ggs)) {
    for (gg in ggs) {
        p <- p + eval(parse(text = gg))
    }
}

png(outfile, width=devpars$width, height=devpars$height, res=devpars$res)
print(p)
dev.off()
