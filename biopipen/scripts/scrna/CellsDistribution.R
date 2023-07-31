source("{{biopipen_dir}}/utils/plot.R")
library(Seurat)
library(rlang)
library(tidyr)
library(dplyr)
library(ggsci)

srtfile = {{in.srtobj | r}}
outdir = {{out.outdir | r}}
{% if in.casefile %}
cases = {{in.casefile | config: "toml" | r: todot="-"}}
{% else %}
cases = {{envs | r: todot="-"}}
{% endif %}

if (length(cases) == 0) {
    stop("No `envs.cases` or `in.casefile` provided.")
}

srtobj = readRDS(srtfile)

mutate_meta = function(meta, mutaters) {
    expr = list()
    for (key in names(mutaters)) {
        expr[[key]] = parse_expr(mutaters[[key]])
    }

    meta %>% mutate(!!!expr)
}

# gg_color_hue <- function(n) {
#     hues = seq(15, 375, length = n + 1)
#     hcl(h = hues, l = 65, c = 100)[1:n]
# }

all_clusters = srtobj@meta.data %>% pull(seurat_clusters)
if (!is.factor(all_clusters)) {
    all_clusters = factor(all_clusters, levels = sort(unique(all_clusters)))
}
# pal = gg_color_hue(nlevels(all_clusters))
# names(pal) = levels(all_clusters)


do_case = function(case) {
    print(paste("- Running for case:", case))
    outfile = file.path(outdir, paste0(case, ".png"))

    casepms = cases$cases[[case]]
    grouppms = casepms$group
    clonepms = casepms$cells

    # mutate meta data
    meta = srtobj@meta.data
    if (!is.null(casepms$mutaters)) {
        meta = mutate_meta(meta, casepms$mutaters)
    }
    if (!is.null(casepms$filter)) {
        meta = meta %>% filter(eval(parse(text=casepms$filter)))
    }

    # meta = meta %>% filter(!is.na(!!sym(grouppms$by)), !is.na(!!sym(clonepms$by)))

    if (!is.null(grouppms$order)) {
        meta[[grouppms$by]] = factor(meta[[grouppms$by]], levels = grouppms$order)
        meta = meta[!is.na(meta[[grouppms$by]]),, drop = FALSE]
    }

    meta = meta |> filter(!is.na(!!sym(clonepms$by)))

    # Sizes
    meta = meta %>%
        add_count(!!sym(clonepms$by), name = "CloneSize") %>%
        add_count(!!sym(clonepms$by), !!sym(grouppms$by), name = "CloneGroupSize") %>%
        add_count(!!sym(clonepms$by), !!sym(grouppms$by), seurat_clusters, name = "CloneGroupClusterSize")
    if (!is.null(clonepms$orderby)) {
        meta = meta %>% arrange(eval(parse(text=clonepms$orderby)))
        order = unique(meta[[clonepms$by]])[1:clonepms$n]
        meta = meta %>% filter(!!sym(clonepms$by) %in% order)
        meta[[clonepms$by]] = factor(meta[[clonepms$by]], levels = order)
    }
    nrows = length(unique(meta[[clonepms$by]]))
    ncols = length(unique(meta[[grouppms$by]]))
    if (is.null(casepms$devpars)) {
        casepms$devpars = list(
            res = 100, width = ncols * 100 + 240, height = max(nrows * 100, 600)
        )
    }

    # plot
    plotGG(
        meta,
        "bar",
        list(
            mapping = aes(
                x = sqrt(CloneGroupSize)/2,
                y = CloneSize,
                width = sqrt(CloneGroupSize),
                fill = seurat_clusters
            ),
            stat = "identity",
            position = "fill"
        ),
        c(
            'geom_col(aes(x=sqrt(CloneGroupSize), y=CloneSize), width=.01, position="fill", color = "#888888")',
            'coord_polar("y", start=0)',
            paste0('facet_grid(vars(', clonepms$by, '), vars(', grouppms$by, '), switch="y")'),
            # 'scale_fill_manual(name = "Cluster", values = pal)',
            # 26-color palette
            'scale_fill_ucscgb(name = "Cluster", limits = levels(all_clusters))',
            'theme_void()',
            'theme(plot.margin = unit(c(1,1,1,1), "cm"))',
            'theme(legend.text=element_text(size=8), legend.title=element_text(size=10))'
        ),
        casepms$devpars,
        outfile
    )
}

sapply(names(cases$cases), do_case)
