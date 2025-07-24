{{ biopipen_dir | joinpaths: "utils", "misc.R" | source_r }}
{{ biopipen_dir | joinpaths: "utils", "plot.R" | source_r }}

library(dplyr)
library(tidyr)
library(ggprism)


immfile <- {{ in.immdata | r }}
outdir <- {{ out.outdir | r }}

subject_key <- {{ envs.subject | r }}
group_key <- {{ envs.group | r }}
sample_order <- {{ envs.order | r }}
diag <- {{ envs.diag | r }}
on <- {{ envs.on | r }}

immdata <- readRDS(immfile)

if (is.null(subject_key) || length(subject_key) == 0) {
    stop("`envs.subject` is required.")
}
if (is.null(group_key)) {
    stop("`envs.group` is required.")
}
if (is.null(sample_order) || length(sample_order) == 0) {
    stop("`envs.order` is required.")
}

# qq plot functions
exponent <- function(x) {
    floor(log10(abs(x)))
}

mantissa <- function(x) {
    mant <- log10(abs(x))
    10^(mant - floor(mant))
}

plot_qq <- function(counts, suf1, suf2) {
    df = data.frame(x = counts[[suf1]], y = counts[[suf2]])
    lim = max(df$x, df$y)
    colnames(df) = c(suf1, suf2)

    ggs = c(
        "theme_prism()",
        paste0("xlim(0, ", lim, ")"),
        paste0("ylim(0, ", lim, ")")
    )
    if (diag) {
        ggs = c(ggs, "geom_abline(intercept = 0, slope = 1, color = 'blue3')")
    }

    plotGG(
        df,
        "point",
        list(aes_string(x=suf1, y=suf2)),
        ggs
    )

}

# qq plots
alpha = 0.375
# https://stats.stackexchange.com/a/77909/292553
subjects = immdata$meta[, subject_key, drop = FALSE] %>% distinct(.keep_all = TRUE)
subjects = subjects[complete.cases(subjects), , drop = FALSE]
for (i in seq_len(nrow(subjects))) {
    # Generate a residency table
    # |    CDR3.aa    | Tumor | Normal |
    # | SEABESRWEFAEF | 0     | 10     |
    # | AWEARWGAWGGGR | 21    | 1      |
    # | GREWFQQWFEWF  | 34    | 0      |
    subject_row <- subjects[i, , drop = FALSE]
    subject <- subject_row %>%
        select(all_of(subject_key)) %>%
        as.character() %>%
        paste(collapse = "-")

    groups <- subject_row %>%
        left_join(immdata$meta, by = subject_key) %>%
        pull(group_key)
    groups <- intersect(sample_order, groups)
    if (length(groups) < 2) {
        warning(paste0("Subject doesn't have enough groups:", subject))
        next
    }

    counts <- list()
    props <- list()
    for (group in groups) {
        sample1 <- subject_row %>%
            left_join(immdata$meta) %>%
            filter(.[[group_key]] == group) %>%
            pull(Sample)
        if (length(sample1) == 0) next
        for (s1 in sample1) {
            counts[[group]] <- c(counts[[group]], immdata$data[[sample1]]$Clones)
            if ("Proportion" %in% names(immdata$data[[sample1]])) {
                props[[group]] <- c(props[[group]], immdata$data[[sample1]]$Proportion)
            }
        }
    }
    lens = sapply(counts, length)
    min_idx = which.min(lens)

    norm_counts = list()
    norm_props = list()
    min_counts = counts[[min_idx]]
    min_props = props[[min_idx]]
    # quantile unequal lengths to the minimal length
    for (i in seq_along(lens)) {
        if (i == min_idx) {
            norm_counts[[groups[i]]] = min_counts
            norm_props[[groups[i]]] = min_props
        } else {
            norm_counts[[groups[i]]] = quantile(
                counts[[i]],
                probs = (rank(min_counts)-alpha)/(length(min_counts)+1-2*alpha)
            )
            if (length(props[[i]]) > 0) {
                norm_props[[groups[i]]] = quantile(
                    props[[i]],
                    probs = (rank(min_props)-alpha)/(length(min_props)+1-2*alpha)
                )
            }
        }
    }

    # all combinations of groups
    for (j in 1:(length(groups) - 1)) {
        for (k in (j + 1):length(groups)) {
            if ("Clones" %in% on) {
                qq_p <- plot_qq(norm_counts, groups[j], groups[k])
                qq_png <- file.path(
                    outdir,
                    paste0("qqplot_", subject, "-", j, ".Clones.png")
                )
                png(qq_png, res = 300, height = 2000, width = 2000)
                print(qq_p)
                dev.off()
            }
            if ("Proportion" %in% on) {
                qq_p <- plot_qq(norm_props, groups[j], groups[k])
                qq_png <- file.path(
                    outdir,
                    paste0("qqplot_", subject, "-", j, ".Proportion.png")
                )
                png(qq_png, res = 300, height = 2000, width = 2000)
                print(qq_p)
                dev.off()
            }
        }
    }
}
