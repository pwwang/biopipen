{{ biopipen_dir | joinpaths: "utils", "misc.R" | source_r }}

library(rlang)
library(ggmanh)

infile <- {{in.infile | r}}
outfile <- {{out.outfile | r}}
chrom_col <- {{envs.chrom_col | r}}
pos_col <- {{envs.pos_col | r}}
pval_col <- {{envs.pval_col | r}}
label_col <- {{envs.label_col | r}}
devpars <- {{envs.devpars | r}}
title <- {{envs.title | r}}
ylabel <- {{envs.ylabel | r}}
rescale <- {{envs.rescale | r}}
rescale_ratio_threshold <- {{envs.rescale_ratio_threshold | r}}
signif <- {{envs.signif | r}}
hicolors <- {{envs.hicolors | r}}
thin_n <- {{envs.thin_n | r}}
thin_bins <- {{envs.thin_bins | r}}
zoom <- {{envs.zoom | r}}
zoom_devpars <- {{envs.zoom_devpars | r}}
chroms <- {{envs.chroms | r}}
args <- {{envs.args | r: todot="-"}}

data <- read.table(infile, header=TRUE, sep="\t", stringsAsFactors=FALSE, check.names = FALSE)

# normalize columns
cnames <- colnames(data)
if (is.numeric(chrom_col)) { chrom_col <- cnames[chrom_col] }
if (is.numeric(pos_col)) { pos_col <- cnames[pos_col] }
if (is.numeric(pval_col)) { pval_col <- cnames[pval_col] }
if (is.numeric(label_col)) { label_col <- cnames[label_col] }

# normalize chroms
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

log_info("Preparing data for plotting ...")
if (length(chroms) == 1 && chroms == "auto") {
    chroms <- unique(data[[chrom_col]])
} else {
    chroms <- norm_chroms(chroms)
}

# prepare data
mp_prep_args = list()
if (length(signif) == 1 && is.character(signif)) {
    signif <- as.numeric(trimws(unlist(strsplit(signif, ","))))
}
siglevel <- min(signif)
if (!is.null(label_col)) {
    data$.label <- ifelse(data[[pval_col]] < siglevel, data[[label_col]], "")
}
if (!is.null(hicolors)) {
    sig_str <- "Significant"
    nsig_str <- "Not significant"
    data$.highlight <- ifelse(data[[pval_col]] < siglevel, sig_str, nsig_str)
    if (length(hicolors) == 1) { hicolors <- c(hicolors, "grey") }
    names(hicolors) <- c(sig_str, nsig_str)
    mp_prep_args$highlight.colname <- ".highlight"
    mp_prep_args$highlight.col <- hicolors
}
mp_prep_args$x <- data
mp_prep_args$chr.colname <- chrom_col
mp_prep_args$pos.colname <- pos_col
mp_prep_args$pval.colname <- pval_col
mp_prep_args$chr.order <- chroms
if (!is.null(thin_n) && thin_n > 0) {
    mp_prep_args$thin.n <- thin_n
    mp_prep_args$thin.bins <- thin_bins
}

mpdata <- do_call(manhattan_data_preprocess, mp_prep_args)

# plot
log_info("Plotting Manhattan plot ...")
args$x <- mpdata
args$signif <- signif
args$plot.title <- title
args$rescale <- rescale
args$rescale.ratio.threshold <- rescale_ratio_threshold
args$y.label <- ylabel
if (!is.null(hicolors)) { args$color.by.highlight <- TRUE }
if (!is.null(label_col)) { args$label.colname <- ".label" }
g <- do_call(manhattan_plot, args)

png(outfile, width=devpars$width, height=devpars$height, res=devpars$res)
print(g)
dev.off()

# zoom into chromosomes
all_chroms <- as.character(unique(mpdata$data[[mpdata$chr.colname]]))
if (!is.null(zoom)) {
    log_info("Zooming into chromosomes ...")
    zoom <- norm_chroms(zoom)
    for (z in zoom) {
        if (!z %in% all_chroms) {
            log_warn("- {z}: not found in data")
            next
        }
        log_info("- {z}")
        args_z <- args
        args_z$chromosome <- z
        args_z$plot.title <- paste0(title, " (", z, ")")
        args_z$x.label <- "Position"
        g_z <- do_call(manhattan_plot, args_z)
        outfile_z <- gsub("\\.png$", paste0("-", z, ".png"), outfile)
        zm_devpars <- zoom_devpars
        zm_devpars$res <- zm_devpars$res %||% devpars$res
        zm_devpars$height <- zm_devpars$height %||% devpars$height
        png(
            outfile_z,
            width=zm_devpars$width,
            height=zm_devpars$height,
            res=zm_devpars$res
        )
        print(g_z)
        dev.off()
    }
}
