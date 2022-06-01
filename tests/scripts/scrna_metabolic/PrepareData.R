library(Seurat)

# Download data (tcell rds and metadata) from
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE123813

meta_file <- {{in.metafile | quote}}
count_file <- {{in.countfile | quote}}
outfile <- {{out.outfile | quote}}
seed = {{envs.seed}}

set.seed(seed)

counts <- read.table(count_file, header = TRUE, row.names = 1, sep = "\t", check.names = F)
metadata <- read.table(meta_file, header = TRUE, row.names = 1, sep = "\t", check.names = F)

# Subset 1000 cells for just testing purpose
counts = counts[, sample(1:ncol(counts), 1000)]
metadata = metadata[colnames(counts),]

# Create seurat object
seurat_obj = CreateSeuratObject(counts=counts)
seurat_obj@meta.data = cbind(
    seurat_obj@meta.data,
    metadata[rownames(seurat_obj@meta.data),]
)
seurat_obj = NormalizeData(seurat_obj)
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all.genes)
seurat_obj <- FindVariableFeatures(object = seurat_obj)
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))

seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

saveRDS(seurat_obj, file=outfile)
