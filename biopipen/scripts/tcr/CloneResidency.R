source("{{biopipen_dir}}/utils/misc.R")

library(rlang)
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
sample_groups = {{ envs.sample_groups | r }}
mutaters = {{ envs.mutaters | r }}
cases = {{ envs.cases | r }}

immdata = readRDS(immfile)
meta = immdata$meta

if (!is.null(mutaters) && length(mutaters) > 0) {
    print("Applying the mutaters ...")
    expr = list()
    for (key in names(mutaters)) {
        expr[[key]] = parse_expr(mutaters[[key]])
    }
    meta = meta %>% mutate(!!!expr)
}

# Fill up cases using `envs.xxx` if not provided and compose a DEFAULT case
# if no cases are provided
if (is.null(cases) || length(cases) == 0) {
    cases = list(
        DEFAULT = list(
            subject = subject_key,
            group = group_key,
            order = sample_order,
            sample_groups = sample_groups
        )
    )
} else {
    for (key in names(cases)) {
        if (is.null(cases[[key]]$subject)) {
            cases[[key]]$subject = subject_key
        }
        if (is.null(cases[[key]]$group)) {
            cases[[key]]$group = group_key
        }
        if (is.null(cases[[key]]$order)) {
            cases[[key]]$order = sample_order
        }
        if (is.null(cases[[key]]$sample_groups)) {
            cases[[key]]$sample_groups = sample_groups
        }
    }
}
# Scatter plot functions
exponent <- function (x) {
    floor(log10(abs(x)))
}

mantissa <- function (x) {
    mant <- log10(abs(x))
    10^(mant - floor(mant))
}

perpare_case = function(casename, case) {
    print(paste("Preparing case:", casename, "..."))
    # Check if required keys are provided
    if (is.null(case$subject) || length(case$subject) == 0) {
        stop(paste("`subject` is required for case:", casename))
    }
    if (is.null(case$group) || length(case$group) == 0) {
        stop(paste("`group` is required for case:", casename))
    }
    if (is.null(case$order) || length(case$order) == 0) {
        stop(paste("`order` is required for case:", casename))
    }

    # Create case-specific directories
    # Scatter plots
    scatter_dir = file.path(outdir, casename, "scatter")
    dir.create(scatter_dir, recursive = TRUE, showWarnings = FALSE)

    # Venn diagrams/Upset plots
    venn_dir = file.path(outdir, casename, "venn")
    dir.create(venn_dir, recursive = TRUE, showWarnings = FALSE)

    # Counts
    counts_dir = file.path(outdir, casename, "counts")
    dir.create(counts_dir, recursive = TRUE, showWarnings = FALSE)
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

handle_subject = function(i, subjects, casename, case) {
    # Generate a residency table
    # |    CDR3.aa    | Tumor | Normal |
    # | SEABESRWEFAEF | 0     | 10     |
    # | AWEARWGAWGGGR | 21    | 1      |
    # | GREWFQQWFEWF  | 34    | 0      |
    subject_row = subjects[i, , drop=FALSE]
    subject = subject_row %>%
        select(all_of(case$subject)) %>%
        as.character() %>%
        paste(collapse="-")

    print(paste("- Handling", i, case$subject, "..."))

    groups = subject_row %>%
        left_join(meta, by=case$subject) %>%
        pull(case$group)
    groups = intersect(case$order, groups)
    if (length(groups) < 2) {
        warning(paste0("-", casename, ", Subject doesn't have enough groups:", subject), immediate. = TRUE)
        return()
    }

    counts = list()
    for (group in groups) {
        sample1 = subject_row %>%
            left_join(meta) %>%
            filter(.[[case$group]] == group) %>%
            pull(Sample)
        if (length(sample1) == 0) next
        if (length(sample1) > 1) {
            warning(paste0(casename, ", Group ", group, " is not unique for subject:", subject), immediate. = TRUE)
        }
        for (s in sample1) {
            counts[[group]] = bind_rows(
                counts[[group]],
                immdata$data[[s]][, c("Clones", "CDR3.aa")]
            )
        }
        counts[[group]] = counts[[group]] %>% mutate(Group = group)
    }
    counts = do_call(bind_rows, counts) %>% pivot_wider(
        id_cols = CDR3.aa,
        names_from = Group,
        values_from = Clones,
        values_fn = function(x) mean(x, na.rm=TRUE)
    )
    counts[is.na(counts)] = 0

    # Save samples to group_by so they can be aligned accordingly in the report
    if (!is.null(sample_groups)) {
        group_dir = file.path(outdir, casename, "sample_groups")
        dir.create(group_dir, showWarnings = FALSE)

        sgroups = subject_row %>%
            left_join(meta) %>%
            pull(sample_groups) %>%
            unique() %>%
            paste(collapse = "-")
        group_file = file.path(group_dir, paste0(sgroups, ".txt"))
        cat(subject, file = group_file, sep = "\n", append = TRUE)
    }

    # Save counts
    counts_dir = file.path(outdir, casename, "counts")
    write.table(
        counts,
        file=file.path(counts_dir, paste0(subject, ".txt")),
        sep="\t",
        row.names=TRUE,
        col.names=TRUE,
        quote=FALSE
    )

    # scatter plot
    # Make plots B ~ A, C ~ B, and C ~ A for order A, B, C
    combns = combn(groups, 2, simplify=FALSE)
    scatter_dir = file.path(outdir, casename, "scatter")
    for (j in seq_along(combns)) {
        pair = combns[[j]]
        scatter_p = plot_scatter(counts, subject, pair[1], pair[2])
        scatter_png = file.path(
            scatter_dir,
            paste0("scatter_", subject, "_", pair[1], "_", pair[2], ".png")
        )
        png(scatter_png, res=300, height=2000, width=2500)
        print(scatter_p)
        dev.off()
    }

    # upset/venn
    venn_dir = file.path(outdir, casename, "venn")
    if (length(groups) > 3) {
        venn_p = counts %>%
            mutate(across(groups), ~ if_else(.x == 0, NA, cur_column())) %>%
            unite(Groups, na.rm=TRUE, sep=":::") %>%
            rowwise() %>%
            mutate(Groups = strsplit(Groups, ":::", fixed = TRUE))
            group_by(CDR3.aa)

        venn_png = file.path(venn_dir, paste0("upset_", subject, ".png"))
    } else {
        venn_data = list()
        for (group in groups) {
            venn_data[[group]] = counts %>% filter(.[[group]] > 0) %>% pull(CDR3.aa)
        }
        venn_p = ggVennDiagram(venn_data)
        venn_png = file.path(venn_dir, paste0("venn_", subject, ".png"))
    }
    png(venn_png, res=300, height=2000, width=2000)
    print(venn_p)
    dev.off()
}

handle_case = function(casename, case) {
    perpare_case(casename, case)
    subjects = meta[, case$subject, drop=F] %>% distinct() %>% drop_na()
    sapply(seq_len(nrow(subjects)), handle_subject, subjects, casename, case)
}

for (casename in names(cases)) {
    handle_case(casename, cases[[casename]])
}
