library(scDblFinder)

# also refer to https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/scDblFinder.html
log_info("Running scDblFinder ...")

if (is.null(envs$scDblFinder$ncores)) {
    envs$scDblFinder$ncores <- envs$ncores
}

envs$scDblFinder$sce <- GetAssayData(sobj, layer = "counts")
if (envs$scDblFinder$ncores > 1) {
    envs$scDblFinder$BPPARAM <- BiocParallel::MulticoreParam(envs$scDblFinder$ncores, RNGseed = 8525)
}
envs$scDblFinder$returnType <- "table"
envs$scDblFinder$ncores <- NULL

doublets <- do_call(scDblFinder, envs$scDblFinder)
doublets <- doublets[doublets$type == "real", , drop = FALSE]
doublets$Barcode <- rownames(doublets)
doublets <- doublets[, c("Barcode", "score", "class"), drop = FALSE]
colnames(doublets) <- c("Barcode", "scDblFinder_score", "scDblFinder_DropletType")
write.table(
    doublets,
    file.path(joboutdir, "scDblFinder_doublets_singlets.txt"),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
)

summary <- as.data.frame(table(doublets$scDblFinder_DropletType))
colnames(summary) <- c("Classification", "Droplet_N")
write.table(
    summary,
    file.path(joboutdir, "scDblFinder_summary.txt"),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
)

log_info("- Adding scDblFinder results to the seurat object ...")
rownames(doublets) <- doublets$Barcode
doublets$Barcode <- NULL
sobj <- AddMetaData(sobj, as.data.frame(doublets))

# Do a dimplot
log_info("- Plotting dimension reduction ...")
dimp <- DimPlot(
    sobj, group.by = "scDblFinder_DropletType", order = "doublet",
    cols = c("#333333", "#FF3333"), pt.size = 0.8, alpha = 0.5)
ggsave(dimp, filename = file.path(plotsdir, "scDblFinder_dimplot.png"))

log_info("- Filtering doublets ...")
sobj <- subset(sobj, cells = doublets$Barcode[doublets$scDblFinder_DropletType == "singlet"])

add_report(
    list(
        kind = "descr",
        content = "The table contains the number of cells classified as singlets and doublets."
    ),
    list(
        kind = "table",
        data = list(path = file.path(joboutdir, "scDblFinder_summary.txt"))
    ),
    h1 = "scDblFinder Results",
    h2 = "The scDblFinder Summary"
)
add_report(
    list(
        name = "Dimension Reduction Plot",
        src = file.path(plotsdir, "scDblFinder_dimplot.png")
    ),
    ui = "table_of_images",
    h1 = "scDblFinder Results",
    h2 = "Plots"
)
