library(dplyr)
library(HGNChelper)
library(Seurat)
library(rlang)

{{ biopipen_dir | joinpaths: "scripts", "scrna", "sctype.R" | source_r }}

sobjfile = {{in.sobjfile | r}}
outfile = {{out.outfile | r}}
tissue = {{envs.sctype_tissue | r}}
db = {{envs.sctype_db | r}}
newcol = {{envs.newcol | r}}
merge_same_labels = {{envs.merge | r}}

if (is.null(db)) { stop("`envs.sctype_args.db` is not set") }

log_info("Reading Seurat object...")
sobj = readRDS(sobjfile)

# prepare gene sets
log_info("Preparing gene sets...")
gs_list = gene_sets_prepare(db, tissue)

scRNAseqData = GetAssayData(sobj, layer = "scale.data")
idents = as.character(unique(Idents(sobj)))
idents = idents[order(as.numeric(idents))]

log_info("Working on different levels of cell type labels ...")
cell_types_list = list()
for (i in seq_along(gs_list)) {
    log_info("- Working on level {i} ...")
    if (is.null(gs_list[[i]])) next

    log_info("  Calculating cell-type scores ...")
    es.max = sctype_score(
        scRNAseqData = scRNAseqData,
        scaled = TRUE,
        gs = gs_list[[i]]$gs_positive,
        gs2 = gs_list[[i]]$gs_negative
    )

    log_info("  Merging cell-type scores by cluster ...")
    cl_resutls = do_call(
        "rbind",
        lapply(
            idents,
            function(cl) {
                es.max.cl = sort(rowSums(es.max[, rownames(sobj@meta.data[sobj@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
                head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(sobj@meta.data$seurat_clusters==cl)), 10)
            }
        )
    )

    sctype_scores = cl_resutls %>%
        group_by(cluster) %>%
        slice_max(scores, n=1, with_ties=TRUE)

    if (nrow(sctype_scores) > length(idents)) {
        sctype_scores_count = sctype_scores %>% count(cluster) %>% filter(n > 1)
        write("\n########## sctype_scores ###########", stderr())
        write(capture.output(sctype_scores), stderr())
        write("\n####### sctype_scores_count ########", stderr())
        write(capture.output(sctype_scores_count), stderr())
        write("\n####################################", stderr())
        log_info("  Scores tied in the above clusters.", immediate. = TRUE)
    }

    if (length(gs_list) == 1 || i > 1) {
        # set low-confident (low ScType score) clusters to "unknown"
        log_info("  Setting low-confident clusters to 'Unknown'...")
        sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"
    }

    celltypes = sapply(
        idents,
        function (cl) {
            cl_type = sctype_scores[sctype_scores$cluster == cl, ]
            as.character(cl_type$type[1])
        }
    )
    names(celltypes) = idents
    cell_types_list[[i]] = celltypes
}

if (length(cell_types_list) == 1) {
    celltypes = cell_types_list[[1]]
} else {
    log_info("Merging cell types at all levels ...")
    celltypes = list()

    for (i in idents) {
        celltypes[[i]] = ""
        for (j in seq_along(cell_types_list)) {
            idt = cell_types_list[[j]][[i]]
            if (idt != "Unknown") {
                celltypes[[i]] = paste(celltypes[[i]], idt)
            }
        }
    }
}


log_info("Renaming cell types...")
ct_numbering = list()
for (key in names(celltypes)) {
    ct = celltypes[[key]]
    ct_numbering[[ct]] = ct_numbering[[ct]] %||% 0
    if (ct_numbering[[ct]] > 0) {
        celltypes[[key]] = paste0(ct, ".", ct_numbering[[ct]])
    }
    ct_numbering[[ct]] = ct_numbering[[ct]] + 1
}

celltypes = as.list(celltypes)
if (is.null(newcol)) {
    sobj$seurat_clusters_id = sobj$seurat_clusters
    celltypes$object = sobj
    sobj = do_call(RenameIdents, celltypes)
    sobj$seurat_clusters = Idents(sobj)
} else {
    celltypes$object = sobj
    sobj = do_call(RenameIdents, celltypes)
    sobj[[newcol]] = Idents(sobj)
    Idents(sobj) = "seurat_clusters"
}
celltypes$object = NULL
gc()

if (merge_same_labels) {
    log_info("Merging clusters with the same labels...")
    sobj <- merge_clusters_with_same_labels(sobj, newcol)
    celltypes <- lapply(celltypes, function(ct) {
        sub("\\.\\d+$", "", ct)
    })
}

log_info("Saving the mappings ...")
write.table(
    data.frame(
        Cluster = names(celltypes),
        Celltype = unlist(celltypes),
        stringsAsFactors = FALSE
    ),
    file = file.path(dirname(outfile), "cluster2celltype.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
)

log_info("Saving Seurat object...")
saveRDS(sobj, outfile)
