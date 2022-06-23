source("{{biopipen_dir}}/utils/misc.R")

library(dplyr)
library(tidyr)
library(ggplot2)
library(ggprism)
library(ggVennDiagram)

theme_set(theme_prism())


immfile = {{ in.immdata | quote }}
outdir = {{ out.outdir | quote }}

subject_key = {{ envs.subject | r }}
group_key = {{ envs.group | r }}
sample_order = {{ envs.order | r }}

immdata = readRDS(immfile)

if (is.null(subject_key) || length(subject_key) == 0) {
    stop("`envs.subject` is required.")
}
if (is.null(group_key)) {
    stop("`envs.group` is required.")
}
if (is.null(sample_order) || length(sample_order) == 0) {
    stop("`envs.order` is required.")
}

# Scatter plot functions
exponent <- function (x) {
    floor(log10(abs(x)))
}

mantissa <- function (x) {
    mant <- log10(abs(x))
    10^(mant - floor(mant))
}

plot_scatter = function(counts, subject, suf1, suf2) {
    # options(repr.plot.width=9, repr.plot.height=7)

    dual = which(counts[[suf1]] > 0 & counts[[suf2]] > 0)
    if (length(dual) <= 2) {
        test = list(estimate=NA, p.value=NA)
    } else {
        test <- cor.test(log(counts[[suf1]][dual]),log(counts[[suf2]][dual]))
    }
    sum_counts1 <- sum(counts[[suf1]])
    sum_counts2 <- sum(counts[[suf2]])

    counts1.norm <- jitter(1+counts[[suf1]], amount=0.25)/sum_counts1
    counts2.norm <- jitter(1+counts[[suf2]], amount=0.25)/sum_counts2

    oo <- sample(length(counts1.norm))
    plotdata = data.frame(x=counts1.norm[oo], y=counts2.norm[oo])
    # plotdata$color = cl.colors[oo]
    names(plotdata) = c(suf1, suf2)
    plotdata = plotdata %>% mutate(
        Type = case_when(
            counts[[suf1]][oo] == 1 & counts[[suf2]][oo] == 0 ~ paste(suf1, "Singleton"),
            counts[[suf1]][oo] == 0 & counts[[suf2]][oo] == 1 ~ paste(suf2, "Singleton"),
            counts[[suf1]][oo] > 1 & counts[[suf2]][oo] == 0 ~ paste(suf1, "Multiplet"),
            counts[[suf1]][oo] == 0 & counts[[suf2]][oo] > 1 ~ paste(suf2, "Multiplet"),
            counts[[suf1]][oo] < counts[[suf2]][oo] ~ "Expanded",
            counts[[suf1]][oo] > counts[[suf2]][oo] ~ "Collapsed",
            TRUE ~ "Dual"
        ),
        Type = as.factor(Type),
        Size = pmax(counts1.norm[oo], counts2.norm[oo])
    )

    xbreaks = c(1/sum_counts1, 0.001+1/sum_counts1, 0.01+1/sum_counts1, 0.1+1/sum_counts1)
    ybreaks = c(1/sum_counts2, 0.001+1/sum_counts2, 0.01+1/sum_counts2, 0.1+1/sum_counts2)

    minx = min(plotdata[[suf1]])
    miny = min(plotdata[[suf2]])
    maxx = max(plotdata[[suf1]])
    maxy = max(plotdata[[suf2]])
    # color = plotdata$color
    # names(color) = color
    # patient = as.character(patient)
    n.formatted <- formatC(length(oo), format="f", big.mark=",", digits=0)
    r.formatted <- format(test$estimate,digits=2,scientific=F)
    if (is.na(test$p.value)) {
        subtitle = bquote(italic(n)[D] == .(length(dual)) ~~ italic(r) == .(r.formatted) ~~ italic(P) == "NA")
    } else if (test$p.value < 1e-4) {
        P.mant <- format(mantissa(test$p.value),digits=2)
        P.exp <- exponent(test$p.value)
        subtitle = bquote(italic(n)[D] == .(length(dual)) ~~ italic(r) == .(r.formatted) ~~ italic(P) == .(P.mant) %*% 10^.(P.exp))
    } else {
        P.formatted <- format(test$p.value,digits=2)
        subtitle = bquote(italic(n)[D] == .(length(dual)) ~~ italic(r) == .(r.formatted) ~~ italic(P) == .(P.formatted))
    }
    ggplot(plotdata) +
        geom_point(
            aes_string(x=bQuote(suf1), y=bQuote(suf2), color="Type", size="Size", fill="Type"),
            alpha=.6,
            shape=21
        ) +
        # geom_point(aes_string(x=x, y=y, color='color'), shape=1) +
        # scale_color_manual(values=color) +
        scale_x_continuous(
            trans="log2",
            limits=c(minx, maxx),
            breaks=xbreaks,
            labels=c("0","0.001","0.01","0.1")
        ) +
        scale_y_continuous(
            trans="log2",
            limits=c(miny, maxy),
            breaks=ybreaks,
            labels=c("0","0.001","0.01","0.1")
        ) +
        theme_prism(base_size = 16) +
        scale_size(guide="none") +
        # theme(legend.position = "none") +
        labs(
            title=bquote(.(subject)~(italic(n) == .(n.formatted))),
            subtitle=subtitle
        ) +
        geom_segment(
            data=data.frame(
                x=c(
                    1.5/sum_counts1,
                    minx,
                    1.5/sum_counts1,
                    minx,
                    2.5/sum_counts1
                ),
                xend=c(
                    maxx, # diagnal
                    maxx,  # horizontal
                    1.5/sum_counts1,   # vertical
                    1.5/sum_counts1,  # horizontal short
                    2.5/sum_counts1  # vertical short

                ),
                y=c(
                    1.5/sum_counts2,
                    1.5/sum_counts2,
                    miny,
                    2.5/sum_counts2,
                    miny
                ),
                yend=c(
                    maxy,
                    1.5/sum_counts2,
                    maxy,
                    2.5/sum_counts2,
                    1.5/sum_counts2
                )
            ),
            aes(x=x, y=y, xend=xend, yend=yend), color='gray'
        )

}

# Scatter plots
scatter_dir = file.path(outdir, "scatter")
dir.create(scatter_dir, showWarnings = FALSE)

# Venn diagrams/Upset plots
venn_dir = file.path(outdir, "venn")
dir.create(venn_dir, showWarnings = FALSE)

# Pie charts
pie_dir = file.path(outdir, "pie")
dir.create(pie_dir, showWarnings = FALSE)

subjects = immdata$meta[, subject_key, drop=F] %>%
    distinct(.[, subject_key], .keep_all = TRUE)
for (i in seq_len(nrow(subjects))) {
    # Generate a residency table
    # |    CDR3.aa    | Tumor | Normal |
    # | SEABESRWEFAEF | 0     | 10     |
    # | AWEARWGAWGGGR | 21    | 1      |
    # | GREWFQQWFEWF  | 34    | 0      |
    subject_row = subjects[i, , drop=FALSE]
    subject = subject_row %>%
        select(all_of(subject_key)) %>%
        as.character() %>%
        paste(collapse="-")

    groups = subject_row %>%
        left_join(immdata$meta, by=subject_key) %>%
        pull(group_key)
    groups = intersect(sample_order, groups)
    if (length(groups) < 2) {
        warning(paste0("Subject doesn't have enough groups:", subject))
        next
    }

    counts = list()
    for (group in groups) {
        sample1 = subject_row %>%
            left_join(immdata$meta) %>%
            filter(.[[group_key]] == group) %>%
            pull(Sample)
        if (length(sample1) == 0) next
        if (length(sample1) > 1) {
            stop(paste0("Group ", suf1, " is not unique for subject:", subject))
        }
        counts[[group]] = as.data.frame(
            immdata$data[[sample1]][, c("Clones", "CDR3.aa")]
        ) %>% mutate(Group = group)
    }
    counts = do_call(bind_rows, counts) %>% pivot_wider(
        CDR3.aa,
        names_from = Group,
        values_from = Clones,
        values_fn = function(x) mean(x, na.rm=TRUE)
    )
    counts[is.na(counts)] = 0

    # scatter plot
    for (j in seq_along(groups)) {
        if (j == 1) next
        scatter_p = plot_scatter(counts, subject, groups[j-1], groups[j])
        scatter_png = file.path(
            scatter_dir,
            paste0("scatter_", subject, ".png")
        )
        png(scatter_png, res=300, height=2000, width=2500)
        print(scatter_p)
        dev.off()
    }

    # upset/venn
    if (length(groups) > 3) {
        venn_p = counts %>%
            mutate(across(groups), ~ if_else(.x == 0, NA, cur_column())) %>%
            unite(Groups, na.rm=TRUE, sep=":::") %>%
            rowwise() %>%
            mutate(Groups = strsplit(Groups, ":::", fixed = TRUE))
            group_by(CDR3.aa)
        venn_png = file.path(
            venn_dir,
            paste0("upset_", subject, ".png")
        )
        # TODO Upset
    } else {

        venn_data = list()
        for (group in groups) {
            venn_data[[group]] = counts %>% filter(.[[group]] > 0) %>% pull(CDR3.aa)
        }
        venn_p = ggVennDiagram(venn_data)
        venn_png = file.path(
            venn_dir,
            paste0("venn_", subject, ".png")
        )
    }
    png(venn_png, res=300, height=2000, width=2000)
    print(venn_p)
    dev.off()

    # Pie chart

}


