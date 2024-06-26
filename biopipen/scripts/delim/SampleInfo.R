source("{{biopipen_dir}}/utils/misc.R")
source("{{biopipen_dir}}/utils/mutate_helpers.R")
library(rlang)
library(dplyr)
library(ggplot2)
library(ggprism)
library(ggrepel)

infile <- {{in.infile | r}}
outfile <- {{out.outfile | r}}
sep <- {{envs.sep | r}}
mutaters <- {{envs.mutaters | r}}
save_mutated <- {{envs.save_mutated | r}}
defaults <- {{envs.defaults | r}}
stats <- {{envs.stats | r}}
exclude_cols <- {{envs.exclude_cols | r}}

if (is.null(exclude_cols)) {
    exclude_cols <- c()
} else {
    exclude_cols <- trimws(unlist(strsplit(exclude_cols, ",")))
}

outdir <- dirname(outfile)
indata <- read.delim(infile, sep = sep, header = TRUE, row.names = NULL)

if (colnames(indata)[1] == "row.names") {
    stop("Wrong number of column names. Do you have the right `sep`?")
}

if (!is.null(mutaters) && length(mutaters) > 0) {
    mutdata <- indata %>%
        mutate(!!!lapply(mutaters, parse_expr))
} else {
    mutdata <- indata
}

write.table(
    if (save_mutated) mutdata else indata,
    file = outfile,
    sep = sep,
    row.names = FALSE,
    col.names = TRUE,
    quote = FALSE
)
add_report(
    list(
        kind = "descr",
        content = "The samples used in the analysis. Each row is a sample, and columns are the meta information about the sample. This is literally the input sample information file, but the paths to the scRNA-seq and scTCR-seq data are hidden.",
        once = TRUE
    ),
    list(
        kind = "table",
        pageSize = 50,
        data = list(file = outfile, sep = sep, excluded = exclude_cols),
        src = FALSE
    ),
    h1 = "Sample Information"
)

theme_set(theme_prism())
for (name in names(stats)) {
    stat <- list_update(defaults, stats[[name]])
    plotfile <- file.path(outdir, paste0(name, ".png"))

    is_continuous <- FALSE
    if (!is.null(stat$subset)) {
        data <- mutdata %>% filter(!!parse_expr(stat$subset))
    } else {
        data <- mutdata
    }
    if (!is.null(stat$group) && !stat$na_group) {
        data <- data %>% filter(!is.na(!!sym(stat$group)))
    }
    if (!is.null(stat$each) && !stat$na_each) {
        data <- data %>% filter(!is.na(!!sym(stat$each)))
    }

    if (is.numeric(data[[stat$on]])) {
        is_continuous <- TRUE
    }

    if (is.null(stat$plot)) {
        stat$plot <- if (is_continuous) "boxplot" else "pie"
    }

    data$..group <- "All"
    group <- if (is.null(stat$group)) sym("..group") else sym(stat$group)
    count_on <- paste0("..count.", stat$on)
    if (!is_continuous) {
        if (!is.null(stat$each)) {
            data <- data %>% add_count(!!group, !!sym(stat$each), name = count_on)
        } else {
            data <- data %>% add_count(!!group, name = count_on)
        }
    }

    if (is.null(stat$devpars)) {
        stat$devpars <- list()
    }
    if (is.null(stat$devpars$width)) {
        stat$devpars$width <- 800
    }
    if (is.null(stat$devpars$height)) {
        stat$devpars$height <- 600
    }
    if (is.null(stat$devpars$res)) {
        stat$devpars$res <- 100
    }

    png(
        plotfile,
        width = stat$devpars$width,
        height = stat$devpars$height,
        res = stat$devpars$res
    )
    if (stat$plot == "boxplot" || stat$plot == "box") {
        p <- ggplot(data, aes(x=!!group, y=!!sym(stat$on), fill=!!group)) +
            geom_boxplot(position = "dodge") +
            scale_fill_biopipen(alpha = .6) +
            xlab("")
    } else if (stat$plot == "violin" ||
               stat$plot == "violinplot" ||
               stat$plot == "vlnplot") {
        p <- ggplot(data, aes(x = !!group, y = !!sym(stat$on), fill=!!group)) +
            geom_violin(position = "dodge") +
            scale_fill_biopipen(alpha = .6) +
            xlab("")
    } else if (
        (grepl("violin", stat$plot) || grepl("vln", stat$plot)) &&
        grepl("box", stat$plot)
    ) {
        p <- ggplot(data, aes(x = !!group, y = !!sym(stat$on), fill = !!group)) +
            geom_violin(position = "dodge") +
            geom_boxplot(width = 0.1, position = position_dodge(0.9), fill="white") +
            scale_fill_biopipen(alpha = .6) +
            xlab("")
    } else if (stat$plot == "histogram" || stat$plot == "hist") {
        p <- ggplot(data, aes(x = !!sym(stat$on), fill = !!group)) +
            geom_histogram(bins = 10, position = "dodge", alpha = 0.8, color = "white") +
            scale_fill_biopipen(alpha = .6)
    } else if (stat$plot == "pie" || stat$plot == "piechart") {
        if (is.null(stat$each)) {
            data <- data %>% distinct(!!group, .keep_all = TRUE)
        } else {
            data <- data %>%
                distinct(!!group, !!sym(stat$each), .keep_all = TRUE) %>%
                mutate(!!group := factor(!!group, levels = unique(!!group))) %>%
                group_by(!!sym(stat$each))
        }
        p <- ggplot(
            data %>% mutate(.size = sum(!!sym(count_on))),
            aes(x = sqrt(.size) / 2, width = sqrt(.size), y = !!sym(count_on), fill = !!group, label = !!sym(count_on))
        ) +
            geom_bar(stat="identity", color="white", position = position_fill(reverse = TRUE)) +
            coord_polar("y", start = 0) +
            theme_void() +
            theme(plot.title = element_text(hjust = 0.5)) +
            geom_label_repel(
                position = position_fill(reverse = TRUE,vjust = .5),
                color="#333333",
                fill="#EEEEEE",
                size=4
            ) +
            scale_fill_biopipen(alpha = .6, name = group) +
            ggtitle(paste0("# ", stat$on))
    } else if (stat$plot == "bar" || stat$plot == "barplot") {
        if (is.null(stat$each)) {
            data <- data %>% distinct(!!group, .keep_all = TRUE)
        } else {
            data <- data %>% distinct(!!group, !!sym(stat$each), .keep_all = TRUE)
        }
        p <- ggplot(
            data,
            aes(x = !!group, y = !!sym(count_on), fill = !!group)) +
            geom_bar(stat = "identity") +
            scale_fill_biopipen(alpha = .6) +
            ylab(paste0("# ", stat$on))
    } else {
        stop("Unknown plot type: ", stat$plot)
    }
    if (!is.null(stat$each)) {
        p <- p + facet_wrap(vars(!!sym(stat$each)), ncol = stat$ncol)
    }
    print(p)
    dev.off()

    by_desc <- ifelse(is.null(stat$by), "", paste0(" by ", stat$by))
    descr <- ifelse(
        is_continuous,
        paste0("The distribution of ", stat$on, by_desc),
        paste0("The number of ", stat$on, by_desc)
    )
    add_report(
        list(kind = "table_image", src = plotfile, name = name, descr = descr),
        h1 = "Statistics",
        ui = "table_of_images:2"
    )
}

save_report(outdir)
