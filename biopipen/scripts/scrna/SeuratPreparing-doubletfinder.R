library(DoubletFinder)

log_info("Running DoubletFinder ...")
log_info("- Preparing Seurat object ...")

if (is.null(envs$DoubletFinder$ncores)) {
    envs$DoubletFinder$ncores <- envs$ncores
}

# More controls from envs?
sobj <- FindNeighbors(sobj, dims = 1:envs$DoubletFinder$PCs)
sobj <- FindClusters(sobj)

log_info("- pK Indentification ...")
sweep.res.list <- paramSweep(
    sobj,
    PCs = 1:envs$DoubletFinder$PCs,
    sct = envs$use_sct,
    num.cores = envs$DoubletFinder$ncores
)
sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
bcmvn <- find.pK(sweep.stats)

bcmvn$Selected <- bcmvn$pK == bcmvn$pK[which.max(bcmvn$BCmetric)[1]]
plot <- ggplot(bcmvn, aes(x = pK, y = BCmetric, color = Selected)) +
    geom_point() +
    # rotate x axis labels
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(plot, filename = file.path(plotsdir, "pK_BCmetric.png"))

pK <- bcmvn$pK[which.max(bcmvn$BCmetric)[1]]
pK <- as.numeric(as.character(pK))
pN <- envs$DoubletFinder$pN
log_info("- Homotypic Doublet Proportion Estimate ...")
homotypic.prop <- modelHomotypic(Idents(sobj))
nExp_poi <- round(nrow(sobj@meta.data) * envs$DoubletFinder$doublets)
nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop))

log_info("- Running DoubletFinder ...")
sobj <- doubletFinder(
    sobj,
    PCs = 1:envs$DoubletFinder$PCs,
    pN = pN,
    pK = pK,
    nExp = nExp_poi.adj,
    reuse.pANN = FALSE,
    sct = envs$use_sct
)
pANN_col <- paste0("pANN_", pN, "_", pK)
pANN_col <- colnames(sobj@meta.data)[grepl(pANN_col, colnames(sobj@meta.data))]
DF_col <- paste0("DF.classifications_", pN, "_", pK)
DF_col <- colnames(sobj@meta.data)[grepl(DF_col, colnames(sobj@meta.data))]
doublets <- as.data.frame(
    cbind(
        colnames(sobj),
        sobj@meta.data[, pANN_col],
        sobj@meta.data[, DF_col]
    )
)
colnames(doublets) <-  c("Barcode","DoubletFinder_score","DoubletFinder_DropletType")
write.table(
    doublets,
    file.path(joboutdir, "DoubletFinder_doublets_singlets.txt"),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
)

summary <- as.data.frame(table(doublets$DoubletFinder_DropletType))
colnames(summary) <- c("Classification", "Droplet_N")
write.table(
    summary,
    file.path(joboutdir, "DoubletFinder_summary.txt"),
    row.names = FALSE,
    quote = FALSE,
    sep = "\t"
)

# Do a dimplot
log_info("- Plotting dimension reduction ...")
dimp <- DimPlot(
    sobj, group.by = DF_col, order = "Doublet",
    cols = c("#333333", "#FF3333"), pt.size = 0.8, alpha = 0.5)
ggsave(dimp, filename = file.path(plotsdir, "DoubletFinder_dimplot.png"))

log_info("- Filtering doublets ...")
sobj <- subset(sobj, cells = doublets$Barcode[doublets$DoubletFinder_DropletType == "Singlet"])

add_report(
    list(
        kind = "descr",
        content = "The table contains the number of cells classified as singlets and doublets."
    ),
    list(
        kind = "table",
        data = list(path = file.path(joboutdir, "DoubletFinder_summary.txt"))
    ),
    h1 = "DoubletFinder Results",
    h2 = "The DoubletFinder Summary"
)
add_report(
    list(
        name = "pK vs BCmetric",
        src = file.path(plotsdir, "pK_BCmetric.png")
    ),
    list(
        name = "Dimension Reduction Plot",
        src = file.path(plotsdir, "DoubletFinder_dimplot.png")
    ),
    ui = "table_of_images",
    h1 = "DoubletFinder Results",
    h2 = "Plots"
)
