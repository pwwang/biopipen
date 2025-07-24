
library(ESCO)
library(rlang)
library(glue)
library(biopipen.utils)

args <- {{envs.esco_args | r: todot="-"}}
args <- args %||% list()

save <- args$save
args$save <- NULL

log <- get_logger()

if (!is.null(seed)) {
    set.seed(seed)
    args$seed <- seed
}
args$nGenes <- ngenes
args$nCells <- nsamples
args$dirname <- paste0(outdir, "/")
args$verbose <- TRUE
args$numCores <- ncores
type <- args$type

log$info("Running simulation ...")
sim <- do_call(escoSimulate, args)
attributes(sim) <- c(attributes(sim), c(simulation_tool = "ESCO"))
save_obj(sim, file.path(outdir, "sim.rds"))

log$info("Plotting ...")
if (type == "single") {
    asys <- assays(sim)
    datalist = list(`simulated-truth` = asys$TrueCounts)
    if (!is.null(asys$counts)) {
        datalist$`zero-inflated` = asys$counts
    }
    if (!is.null(asys$observedcounts)) {
        datalist$`down-sampled` = asys$observedcounts
    }

    log$info("- Plotting the data ...")
    dataplot <- file.path(outdir, "data.png")
    png(dataplot, width=length(datalist) * 600, height=1200, res=30)
    heatdata(datalist, norm = FALSE, size = 2, ncol = 3)
    dev.off()

    rholist <- metadata(sim)$Params@corr
    if (length(rholist) > 0) {
        log$info("- Plotting the GCN ...")
        corrgenes <- rownames(rholist[[1]])
        gcnlist = lapply(datalist, function(data)gcn(data, genes = corrgenes))
        gcnlist = append(gcnlist, list("given truth" = rholist[[1]]), 1)

        gcnplot <- file.path(outdir, "gcn.png")
        png(gcnplot, width=length(gcnlist) * 600, height=1200, res=30)
        heatgcn(gcnlist, size = 2, ncol = 4)
        dev.off()
    }
} else if (type == "groups") {
    asys <- assays(sim)
    # organize the marker gene info
    genegroup = paste0("Group", rowData(sim)$GeneGroup)
    genegroup[which(genegroup=="Group0")] = "None"
    geneinfo = data.frame(genes = rowData(sim)$Gene,
                          newcelltype = as.factor(genegroup))

    # organize the cell info
    cellinfo = data.frame(çells = colData(sim)$Cell,
                          newcelltype= as.factor(colData(sim)$Group))

    # data
    datalist = list(`simulated-truth` = asys$TrueCounts)
    if (!is.null(asys$counts)) {
        datalist$`zero-inflated` = asys$counts
    }
    if (!is.null(asys$observedcounts)) {
        datalist$`down-sampled` = asys$observedcounts
    }

    log$info("- Plotting the data ...")
    dataplot <- file.path(outdir, "data.png")
    png(dataplot, width=length(datalist) * 600, height=1200, res=30)
    heatdata(datalist, cellinfo = cellinfo, geneinfo = geneinfo, size = 1, ncol = 3)
    dev.off()

    log$info("- Plotting the GCN for all marker genes (i.e. DE genes) across all cell groups ...")
    degeneinfo = geneinfo[which(geneinfo$newcelltype!="None"),]
    degeneinfo$newcelltype = droplevels(degeneinfo$newcelltype)
    degcnlist = lapply(datalist, function(data)gcn(data, genes = degeneinfo$genes))
    gcnplot <- file.path(outdir, "gcn-allgroups.png")
    png(gcnplot, width=length(degcnlist) * 700, height=1200, res=30)
    heatgcn(degcnlist, geneinfo = degeneinfo, size = 2, ncol = 3)
    dev.off()

    log$info("- Plotting the GCN for marker genes within one cell group ...")
    rholist = metadata(sim)$Params@corr
    group2_gcnlist = lapply(datalist,
                            function(data){
                            gcn(data[,which(colData(sim)$Group=="Group2")],
                                CPM2 = TRUE,
                                genes = rownames(rholist[["Group2"]]))})
    group2_gcnlist = append(group2_gcnlist,
                            list("given truth" = rholist[["Group2"]]), 1)
    gcnplot2 <- file.path(outdir, "gcn-onegroup.png")
    png(gcnplot2, width=length(group2_gcnlist) * 700, height=1200, res=30)
    heatgcn(group2_gcnlist, size = 3, ncol = 4)
    dev.off()
} else if (type == "tree") {
    # get the data
    datatrue = assays(sim)$TrueCounts

    # get the cellinfo
    cellinfo  = data.frame(cell = colData(sim)$Cell,
                        newcelltype = as.factor(colData(sim)$Group))
    levels(cellinfo$newcelltype) = tree$tip.label

    # get the geneinfo
    genegroup = paste0("Group", rowData(sim)$GeneGroup)
    genegroup[which(genegroup=="Group0")] = "None"
    geneinfo = data.frame(genes = rowData(sim)$Gene,
                        newcelltype = as.factor(genegroup))
    levels(geneinfo$newcelltype)[1:3] = tree$tip.label

    # get the DE geneinfo
    groups <- colData(sim)$Group
    group.names <- sort(unique(groups))
    group.facs.gene <- rowData(sim)[, paste0("DEFac", group.names)]
    DEgene.name = as.character(rowData(sim)$Gene[which(group.facs.gene[,1]>1)])
    degeneinfo = geneinfo[match(DEgene.name, geneinfo$genes),]

    log$info("- Plotting the data ...")
    dataplot <- file.path(outdir, "data.png")
    png(dataplot, width=2000, height=1200, res=30)
    # plot the data
    heatdata(list(datatrue),
            colv = TRUE,
            cellinfo = cellinfo,
            geneinfo = degeneinfo,
            genes = degeneinfo$genes,
            size = 1.5, ncol = 1)
    dev.off()
} else if (type == "traj") {
    datatrue = assays(sim)$TrueCounts

    # get the cellinfo
    cellinfo = data.frame(cell = colData(sim)$Cell,
                        newcelltype = colData(sim)$Path)
    # get the pesudo time
    celltime = data.frame(path = as.numeric(colData(sim)$Path),
                        step = as.numeric(colData(sim)$Step))
    celltime = order(celltime[,1], celltime[,2])

    # get the geneinfo
    degenes = which(metadata(sim)$Params@paths.DEgenes==1)

    log$info("- Plotting the trajectory ...")
    trajplot <- file.path(outdir, "traj.png")
    png(trajplot, width=1600, height=1200, res=30)
    # plot the data
    umapplot(t(t(datatrue)/colSums(datatrue)),
            celltype  = colData(sim)$Path,
            labels = levels(as.factor(colData(sim)$Path)))
    dev.off()

    log$info("- Plotting the data ...")
    dataplot <- file.path(outdir, "data.png")
    heatdata(list("simulated truth" = datatrue[degenes,]),
         cellinfo = cellinfo,
         colv = celltime, size = 1, ncol = 1)
    dev.off()
}

simulated <- switch(save,
    `simulated-truth` = assays(sim)$TrueCounts,
    `zero-inflated` = assays(sim)$counts,
    `down-sampled` = assays(sim)$observedcounts,
    { stop(glue("Unknown save option: {save}, expected one of 'simulated-truth', 'zero-inflated', 'down-sampled'")) }
)

