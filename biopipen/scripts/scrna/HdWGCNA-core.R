############## Setup ##############
log$info("Setting up the WGCNA analysis ...")
srtobj <- do_call("SetupForWGCNA", c(list(seurat_obj = srtobj, wgcna_name = wgcna_name), SetupForWGCNAArgs))

############## Metacells or pseudobulk ##############
if (use_pseudobulk) {
    log$info("Aggregating the cells into pseudobulk ...")
    counts <- Seurat::GetAssayData(srtobj, layer = "counts")
    se <- do_call("AggregatePseudobulk", c(
        list(X = counts, meta = srtobj@meta.data),
        AggregatePseudobulkArgs
    ))
    se <- do_call("NormalizeCounts", c(list(se = se), NormalizeCountsArgs))
} else {
    log$info("Constructing metacells ...")
    srtobj <- do_call("MetacellsByGroups", c(list(seurat_obj = srtobj), MetacellsByGroupsArgs))
    srtobj <- do_call("NormalizeMetacells", c(list(seurat_obj = srtobj), NormalizeMetacellsArgs))
    if (!is.null(ScaleMetacellsArgs)) {
        srtobj <- do_call("ScaleMetacells", c(list(seurat_obj = srtobj), ScaleMetacellsArgs))
    }
    if (!is.null(RunPCAMetacellsArgs)) {
        srtobj <- do_call("RunPCAMetacells", c(list(seurat_obj = srtobj), RunPCAMetacellsArgs))
    }
    if (!is.null(RunHarmonyMetacellsArgs)) {
        if (is.null(RunHarmonyMetacellsArgs$group.by.vars)) {
            stop("`RunHarmonyMetacells` requires `group-by-vars` to be set!")
        }
        srtobj <- do_call("RunHarmonyMetacells", c(list(seurat_obj = srtobj), RunHarmonyMetacellsArgs))
    }
    if (!is.null(RunUMAPMetacellsArgs)) {
        srtobj <- do_call("RunUMAPMetacells", c(list(seurat_obj = srtobj), RunUMAPMetacellsArgs))
    }
}

############## SetDatExpr / SetMultiExpr ##############
if (use_consensus) {
    if (is.null(SetMultiExprArgs)) {
        stop("`SetMultiExpr` must be set when `use_consensus` is true!")
    }
    if (use_pseudobulk) {
        SetMultiExprArgs$mat <- se
        if (is.null(SetMultiExprArgs$layer)) SetMultiExprArgs$layer <- "VST"
    }
    log$info("Setting up the multi-group expression data for the consensus network ...")
    srtobj <- do_call("SetMultiExpr", c(list(seurat_obj = srtobj, wgcna_name = wgcna_name), SetMultiExprArgs))
    log$info("Testing the soft powers for the consensus network ...")
    srtobj <- do_call("TestSoftPowersConsensus", c(list(seurat_obj = srtobj, wgcna_name = wgcna_name), TestSoftPowersConsensusArgs))
    ConstructNetworkArgs$consensus <- TRUE
} else {
    if (use_pseudobulk) {
        SetDatExprArgs$mat <- se
        if (is.null(SetDatExprArgs$layer)) SetDatExprArgs$layer <- "VST"
        SetDatExprArgs$use_metacells <- NULL
        SetDatExprArgs$group_name <- NULL
        SetDatExprArgs$group.by <- NULL
    }
    log$info("Setting up the expression data ...")
    srtobj <- do_call("SetDatExpr", c(list(seurat_obj = srtobj, wgcna_name = wgcna_name), SetDatExprArgs))
    log$info("Testing the soft powers ...")
    srtobj <- do_call("TestSoftPowers", c(list(seurat_obj = srtobj, wgcna_name = wgcna_name), TestSoftPowersArgs))
}

############## ConstructNetwork ##############
tom_dir <- ConstructNetworkArgs$tom_outdir %||% file.path(cache %||% outdir, "TOM")
dir.create(tom_dir, recursive = TRUE, showWarnings = FALSE)
ConstructNetworkArgs$tom_outdir <- tom_dir
if (is.null(ConstructNetworkArgs$tom_name)) {
    ConstructNetworkArgs$tom_name <- wgcna_name
}
log$info(glue("Constructing the network, TOM matrices at {tom_dir} ..."))
srtobj <- do_call("ConstructNetwork", c(list(seurat_obj = srtobj, wgcna_name = wgcna_name), ConstructNetworkArgs))
# hdWGCNA stores the TOM path as paste0(getwd(), "/", tom_outdir), which
# mangles absolute tom_outdir paths; undo the prefix when it does not exist
net <- srtobj@misc[[wgcna_name]]$wgcna_net
tomfiles <- net$TOMFiles
if (length(tomfiles) > 0 && !file.exists(tomfiles[[1]])) {
    fixed <- sub(paste0("^", getwd(), "/"), "", tomfiles)
    log$warn(glue("Fixing the TOM file path: {tomfiles[[1]]} -> {fixed[[1]]}"))
    net$TOMFiles <- fixed
    srtobj@misc[[wgcna_name]]$wgcna_net <- net
}

############## Module eigengenes and connectivity ##############
log$info("Computing the module eigengenes ...")
srtobj <- do_call("ModuleEigengenes", c(list(seurat_obj = srtobj, wgcna_name = wgcna_name), ModuleEigengenesArgs))

log$info("Computing the module connectivity ...")
srtobj <- do_call("ModuleConnectivity", c(list(seurat_obj = srtobj, wgcna_name = wgcna_name), ModuleConnectivityArgs))

if (!is.null(ResetModuleNamesArgs)) {
    log$info("Resetting the module names ...")
    srtobj <- do_call("ResetModuleNames", c(list(seurat_obj = srtobj, wgcna_name = wgcna_name), ResetModuleNamesArgs))
}
if (!is.null(ResetModuleColorsArgs)) {
    log$info("Resetting the module colors ...")
    srtobj <- do_call("ResetModuleColors", c(list(seurat_obj = srtobj, wgcna_name = wgcna_name), ResetModuleColorsArgs))
}

srtobj <- SetActiveWGCNA(srtobj, wgcna_name)

if (!is.null(ModuleExprScoreArgs)) {
    log$info("Computing the module expression scores ...")
    srtobj <- do_call("ModuleExprScore", c(list(seurat_obj = srtobj, wgcna_name = wgcna_name), ModuleExprScoreArgs))
}
