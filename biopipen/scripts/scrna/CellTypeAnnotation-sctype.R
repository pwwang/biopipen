library(dplyr)
library(HGNChelper)
library(Seurat)
library(rlang)
library(biopipen.utils)

{% include biopipen_dir + "/scripts/scrna/sctype.R" %}

sobjfile = {{in.sobjfile | r}}
outfile = {{out.outfile | r}}
tissue = {{envs.sctype_tissue | r}}
db = {{envs.sctype_db | r}}
newcol = {{envs.newcol | r}}
ident = {{envs.ident | r }}
merge_same_labels = {{envs.merge | r}}

if (is.null(db)) { stop("`envs.sctype_args.db` is not set") }

log <- get_logger()

log$info("Reading Seurat object...")
sobj = biopipen.utils::read_obj(sobjfile)
ident <- ident %||% biopipen.utils::GetIdentityColumn(sobj)
Idents(sobj) <- ident

# prepare gene sets
log$info("Preparing gene sets...")
gs_list = gene_sets_prepare(db, tissue)

scRNAseqData = GetAssayData(sobj, layer = "scale.data")
idents = as.character(unique(Idents(sobj)))
idents = idents[order(as.numeric(idents))]

log$info("Working on different levels of cell type labels ...")
cell_types_list = list()
for (i in seq_along(gs_list)) {
    log$info("- Working on level {i} ...")
    if (is.null(gs_list[[i]])) next

    log$info("  Calculating cell-type scores ...")
    es.max = sctype_score(
        scRNAseqData = scRNAseqData,
        scaled = TRUE,
        gs = gs_list[[i]]$gs_positive,
        gs2 = gs_list[[i]]$gs_negative
    )

    log$info("  Merging cell-type scores by cluster ...")
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
        log$info("  Scores tied in the above clusters.", immediate. = TRUE)
    }

    if (length(gs_list) == 1 || i > 1) {
        # set low-confident (low ScType score) clusters to "unknown"
        log$info("  Setting low-confident clusters to 'Unknown'...")
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
    log$info("Merging cell types at all levels ...")
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


log$info("Renaming cell types...")
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
    sobj <- rename_idents(sobj, ident, celltypes)
} else {
    sobj@meta.data[[newcol]] = celltypes[as.character(Idents(sobj))]
}
celltypes$object = NULL
gc()

if (merge_same_labels) {
    log$info("Merging clusters with the same labels...")
    sobj <- merge_clusters_with_same_labels(sobj, newcol)
    celltypes <- lapply(celltypes, function(ct) {
        sub("\\.\\d+$", "", ct)
    })
}

log$info("Saving the mappings ...")
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

log$info("Saving Seurat object...")
biopipen.utils::save_obj(sobj, outfile)
