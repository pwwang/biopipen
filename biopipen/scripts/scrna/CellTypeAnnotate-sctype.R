library(dplyr)
library(HGNChelper)
library(Seurat)
source("{{biopipen_dir}}/scripts/scrna/sctype.R")

sobjfile = {{in.sobjfile | r}}
outfile = {{out.outfile | r}}
tissue = {{envs.sctype_tissue | r}}
db = {{envs.sctype_db | r}}

if (is.null(tissue)) { stop("`envs.sctype_args.tissue` is not set") }
if (is.null(db)) { stop("`envs.sctype_args.db` is not set") }

sobj = readRDS(sobjfile)

# prepare gene sets
gs_list = gene_sets_prepare(db, tissue)

# get cell-type by cell matrix
es.max = sctype_score(
    scRNAseqData = GetAssayData(sobj, slot = "scale.data"),
    scaled = TRUE,
    gs = gs_list$gs_positive,
    gs2 = gs_list$gs_negative
)

# NOTE: scRNAseqData parameter should correspond to your input scRNA-seq matrix.
# In case Seurat is used, it is either sobj[["RNA"]]@scale.data (default), sobj[["SCT"]]@scale.data, in case sctransform is used for normalization,
# or sobj[["integrated"]]@scale.data, in case a joint analysis of multiple single-cell datasets is performed.

idents = as.character(unique(Idents(sobj)))
idents = idents[order(as.numeric(idents))]
# merge by cluster
cl_resutls = do.call(
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
    warning("Scores tied in the above clusters.")
}

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"

celltypes = sapply(
    idents,
    function (cl) {
        cl_type = sctype_scores[sctype_scores$cluster == cl, ]
        as.character(cl_type$type[1])
    }
)
celltypes = make.unique(celltypes)
names(celltypes) = idents
celltypes$object = sobj

sobj = do.call(RenameIdents, celltypes)
sobj$seurat_clusters = Idents(sobj)

saveRDS(sobj, outfile)
