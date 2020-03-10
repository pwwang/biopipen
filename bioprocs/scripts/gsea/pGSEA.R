
{{"__init__.R" | rimport}}
library(fgsea)
library(data.table)

exprfile <- {{i.exprfile | R}}
sifile <- {{i.saminfo | R}}
gmtfile <- {{i.gmtfile | R}}
ranking <- {{args.ranking | R}}
nthread <- {{ args.nthread | R}}
seed <- {{ args.seed | R}}
inopts <- {{ args.inopts | R }}
devpars <- {{args.devpars | R}}
minsize <- {{args.minsize | R}}
top <- {{args.top | R}}
outfile <- {{ o.outfile | R }}
outdir <- {{ o.outdir | R }}

if (seed > 0) {
    set.seed(seed)
}

inopts$cnames <- TRUE
inopts$rnames <- TRUE
exprs <- read.table.inopts(exprfile, inopts)
saminfo <- utils$sampleinfo$SampleInfo2(sifile)
samples <- colnames(exprs)

get_group_samples <- function() {
    groups <- saminfo$all_groups()
    if (length(groups) != 2) {
        stop("Can only handle 2 groups currently.")
    }
    ret <- list()
    ret[[groups[1]]] <- match(saminfo$get_samples(by = "Group",
                                                  value = groups[1]), samples)
    ret[[groups[2]]] <- match(saminfo$get_samples(by = "Group",
                                                  value = groups[2]), samples)
    ret
}

s2n <- function(row, index1, index2) {
    (mean(row[index1], na.rm = TRUE) - mean(row[index2], na.rm = TRUE)) /
    (sd(row[index1], na.rm = TRUE) + sd(row[index2], na.rm = TRUE))
}

ttest <- function(row, index1, index2) {
    (mean(row[index1], na.rm = TRUE) - mean(row[index2], na.rm = TRUE)) /
    sqrt(sd(row[index1], na.rm = TRUE)^2 / length(index1) +
         sd(row[index2], na.rm = TRUE)^2 / length(index2))
}

roc <- function(row, index1, index) {
    mean(row[index1], na.rm = TRUE) / mean(row[index2], na.rm = TRUE)
}

doc <- function(row, index1, index2) {
    mean(row[index1], na.rm = TRUE) - mean(row[index2], na.rm = TRUE)
}

log2roc <- function(row, index1, index2) {
    log2(mean(row[index1], na.rm = TRUE)) -
    log2(mean(row[index2], na.rm = TRUE))
}

gene_ranking <- function() {
    group_samples <- get_group_samples()
    index1 <- group_samples[[1]]
    index2 <- group_samples[[2]]
    ret <- apply(exprs, 1, function(row) {
        # add perturbation to avoid sd = 0
        do.call(ranking, list(row + rnorm(length(row), 0, 1e-6),
                              index1, index2))
    })
    names(ret) <- rownames(exprs)
    sort(ret)
}

mean_center <- function(indata) {
    rmeans <- rowMeans(indata, na.rm = TRUE)
    rstdev <- apply(indata, 1, sd, na.rm = TRUE)
    rstdev[rstdev == 0] <- 1e-6
    2 * (indata - rmeans) / rstdev
}

pathways <- gmtPathways(gmtfile)
ranks <- gene_ranking()
fgsea_res <- fgseaMultilevel(pathways = pathways,
                             stats = ranks,
                             minSize = minsize,
                             maxSize = max(sapply(pathways, length)),
                             #nperm = nperm,
                             nproc = nthread)

setorder(fgsea_res, "pval")
fwrite(fgsea_res, file = outfile, sep = "\t", sep2 = c("", ",", ""))

# plot
if (top > 0) {

    tableplot <- file.path(outdir,
                           paste0(tools::file_path_sans_ext(basename(outfile)),
                                  ".gseatable.png"))
    do.call(png, c(list(tableplot), devpars))
    plotGseaTable(pathways[fgsea_res[1:top, pathway]], ranks,
                  fgsea_res, gseaParam = 0.5)
    dev.off()

    {{"plot.R" | rimport}}
    group_samples <- get_group_samples()
    groups = names(group_samples)
    index1 <- group_samples[[1]]
    index2 <- group_samples[[2]]
    sampleanno <- data.frame(Sample = c(samples[index1], samples[index2]),
                             Group = c(rep(groups[1], length(index1)),
                                       rep(groups[2], length(index2))))
    group_colors <- c()
    group_colors[groups[1]] <- "red3"
    group_colors[groups[2]] <- "green3"
    plot_one <- function(i) {
        pathway <- fgsea_res[i, pathway]
        # enrichment plot
        enplot <- file.path(outdir, paste0(pathway, ".enrich.png"))
        do.call(png, c(list(enplot), list(res = devpars$res,
                                          width = devpars$width,
                                          height = devpars$height / 2)))
        print(plotEnrichment(pathways[[pathway]], ranks, gseaParam = 0.5) +
              ggplot2::labs(title = pathway))
        dev.off()

        # heatmap for leading
        leading <- unlist(fgsea_res[i, leadingEdge])
        plotdata <- exprs[leading, c(index1, index2)]
        hmfile <- file.path(outdir, paste0(pathway, ".leading_heatmap.png"))
        plot.heatmap2(mean_center(plotdata),
                      hmfile,
                      params = list(
                          row_dend_reorder = TRUE,
                          column_title = paste(
                              "Heatmap for leading edge genes in pathway:",
                              pathway
                          ),
                          cluster_rows = length(leading) >= 2,
                          top_annotation = HeatmapAnnotation(
                              Group = sampleanno$Group,
                              col = list(Group = group_colors)
                          ),
                          heatmap_legend_param = list(title = "Expression")
                      ),
                      draw = list(merge_legends = TRUE))

        # heatmap for all genes in pathway
        genes = intersect(pathways[[pathway]], rownames(exprs))
        plotdata <- exprs[genes, , drop = FALSE]
        hmfile <- file.path(outdir, paste0(pathway, ".heatmap.png"))
        plot.heatmap2(mean_center(plotdata),
                      hmfile,
                      params = list(
                          row_dend_reorder = TRUE,
                          show_row_names = nrow(plotdata) <= 50,
                          column_title_gp = gpar(fontsize = 10),
                          column_title = paste(
                              "Heatmap for genes in pathway:",
                              pathway
                          ),
                          row_names_gp = gpar(fontsize = 6),
                          cluster_rows = length(leading) >= 2,
                          top_annotation = HeatmapAnnotation(
                              Group = sampleanno$Group,
                              col = list(Group = group_colors)
                          ),
                          heatmap_legend_param = list(title = "Expression")
                      ),
                      draw = list(merge_legends = TRUE))
    }

    library(parallel)
    iferror <- mcmapply(plot_one, 1:top, mc.cores = nthread)
    if (any(grepl("Error", iferror, fixed = TRUE))) {
        stop(iferror)
    }
}
