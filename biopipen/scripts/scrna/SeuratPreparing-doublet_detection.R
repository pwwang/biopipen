.get_envs_cached_doubletfinder <- function() {
    envs_cache <- envs
    envs_cache$ncores <- NULL
    envs_cache$doublet_detector <- NULL
    envs_cache$scDblFinder <- NULL
    envs_cache$DoubletFinder$ncores <- NULL
    envs_cache
}

.get_envs_cached_scdblfinder <- function() {
    envs_cache <- envs
    envs_cache$ncores <- NULL
    envs_cache$doublet_detector <- NULL
    envs_cache$DoubletFinder <- NULL
    envs_cache$scDblFinder$ncores <- NULL
    envs_cache
}

.run_doubletfinder <- function() {
    library(DoubletFinder)
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
    doublets <- sobj@meta.data[, c(pANN_col, DF_col), drop = FALSE]
    colnames(doublets) <-  c("DoubletFinder_score","DoubletFinder_DropletType")
    doublets$DoubletFinder_DropletType <- tolower(doublets$DoubletFinder_DropletType)

    pk_plot <- ggplot(bcmvn, aes(x = pK, y = BCmetric, color = Selected)) +
        geom_point() +
        # rotate x axis labels
        theme(axis.text.x = element_text(angle = 90, hjust = 1))
    list(doublets = doublets, pk_plot = pk_plot)
}

.run_scdblfinder <- function() {
    library(scDblFinder)
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
    doublets <- doublets[, c("score", "class"), drop = FALSE]
    colnames(doublets) <- c("scDblFinder_score", "scDblFinder_DropletType")

    list(doublets = doublets)
}

run_dd <- function(detector) {
    log_info("Running {detector} ...")
    if (detector == "DoubletFinder") {
        envs_cache_fun <- .get_envs_cached_doubletfinder
        run_fun <- .run_doubletfinder
    } else if (detector == "scDblFinder") {
        envs_cache_fun <- .get_envs_cached_scdblfinder
        run_fun <- .run_scdblfinder
    } else {
        stop("Unknown doublet detector: ", detector)
    }

    cached <- get_cached(envs_cache_fun(), detector, cache_dir)
    if (!is.null(cached$data)) {
        log_info("- Loading cached results ...")
        results <- cached$data
    } else {
        results <- run_fun()

        cached$data <- results
        save_to_cache(cached, detector, cache_dir)
    }

    results
}

save_dd <- function(dd, detector) {
    doublets <- dd$doublets
    write.table(
        doublets,
        file.path(joboutdir, paste0(detector, "_doublets_singlets.txt")),
        row.names = FALSE,
        quote = FALSE,
        sep = "\t"
    )

    summary <- as.data.frame(table(dd$doublets[[paste0(detector, "_DropletType")]]))
    colnames(summary) <- c("Classification", "Droplet_N")
    write.table(
        summary,
        file.path(joboutdir, paste0(detector, "_summary.txt")),
        row.names = FALSE,
        quote = FALSE,
        sep = "\t"
    )

    n_doublet <- summary$Droplet_N[summary$Classification == 'doublet']
    log_info("- {n_doublet}/{sum(summary$Droplet_N)} doublets detected.")
}

add_dd_to_seurat <- function(sobj, dd) {
    AddMetaData(sobj, metadata = as.data.frame(dd$doublets))
}

plot_dd <- function(sobj, dd, detector) {
    if (detector == "DoubletFinder") {
        log_debug("- Plotting pK vs BCmetric ...")
        ggsave(dd$pk_plot, filename = file.path(plotsdir, "DoubletFinder_pK_BCmetric.png"))
    }

    log_info("- Plotting dimension reduction ...")
    dimp <- DimPlot(
        sobj, group.by = paste0(detector, "_DropletType"), order = "doublet",
        cols = c("#333333", "#FF3333"), pt.size = 0.8, alpha = 0.5)
    ggsave(dimp, filename = file.path(plotsdir, paste0(detector, "_dimplot.png")))
}

filter_dd <- function(sobj, dd, detector) {
    subset(sobj,
        cells = rownames(dd$doublets[
            dd$doublets[[paste0(detector, "_DropletType")]] == "singlet", ,
            drop = FALSE
        ]))
}

report_dd <- function(detector) {
    add_report(
        list(
            kind = "descr",
            content = "The table contains the number of cells classified as singlets and doublets."
        ),
        list(
            kind = "table",
            data = list(path = file.path(joboutdir, paste0(detector, "_summary.txt")))
        ),
        h1 = paste0(detector, " Results"),
        h2 = paste0("The ", detector, " Summary")
    )

    if (detector == "DoubletFinder") {
        add_report(
            list(name = "pK vs BCmetric", src = file.path(plotsdir, "pK_BCmetric.png")),
            list(name = "Dimension Reduction Plot", src = file.path(plotsdir, "DoubletFinder_dimplot.png")),
            ui = "table_of_images",
            h1 = "DoubletFinder Results",
            h2 = "Plots"
        )
    } else {
        add_report(
            list(name = "Dimension Reduction Plot",src = file.path(plotsdir, "scDblFinder_dimplot.png")),
            ui = "table_of_images",
            h1 = "scDblFinder Results",
            h2 = "Plots"
        )
    }
}
