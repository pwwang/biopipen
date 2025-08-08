library(glue)
library(dplyr)
library(tidyr)
library(tibble)
library(Seurat)
library(biopipen.utils)

screpdata <- {{in.screpdata | r}}
outfile <- {{out.outfile | r}}
joboutdir <- {{job.outdir | r}}
python <- {{envs.python | r}}
within_sample <- {{envs.within_sample | r}}
assay <- {{envs.assay | r}}
predefined_b <- {{envs.predefined_b | r}}
max_iter <- {{envs.max_iter | int}}
save_tessa <- {{envs.save_tessa | r}}

log <- get_logger()
reporter <- get_reporter()

# In case this script is running in the cloud and <biopipen_dir> can not be found in there
# In stead, we use the python command, which is associated with the cloud environment,
# to get the biopipen directory
biopipen_dir <- get_biopipen_dir(python)
tessa_srcdir <- file.path(biopipen_dir, "scripts", "tcr", "TESSA_source")

outdir <- dirname(outfile)
result_dir <- file.path(outdir, "result")
tessa_dir <- file.path(outdir, "tessa")
if (!dir.exists(result_dir)) dir.create(result_dir)
if (!dir.exists(tessa_dir)) dir.create(tessa_dir)

### Start preparing input files for TESSA
# Prepare input files
log$info("Reading input file ...")
sobj <- read_obj(screpdata)

log$info("Preparing TCR input file ...")
# If immfile endswith .rds, then it is an immunarch object
tcrdata <- sobj@meta.data %>%
    rownames_to_column("contig_id") %>%
    select(contig_id, CTaa, CTgene, sample = Sample) %>%
    filter(!is.na(CTaa) & !is.na(CTgene)) %>%
    separate(CTaa, into = c(NA, "cdr3"), sep = "_", remove = TRUE) %>%
    filter(!is.na(cdr3) & cdr3 != "NA" & cdr3 != "nan") %>%
    separate(CTgene, into = c(NA, "vjgene"), sep = "_", remove = TRUE) %>%
    separate(vjgene, into = c("v_gene", NA, "j_gene", NA), sep = "\\.", remove = TRUE) %>%
    mutate(v_gene = sub("-\\d+$", "", v_gene), j_gene = sub("-\\d+$", "", j_gene))

log$info("Preparing expression input file ...")
expr <- GetAssayData(sobj, layer = "data")
cell_ids <- intersect(tcrdata$contig_id, colnames(expr))
# Warning about unused cells
unused_expr_cells <- setdiff(colnames(expr), cell_ids)
if (length(unused_expr_cells) > 0) {
    log$warn(glue("{length(unused_expr_cells)}/{ncol(expr)} cells without TCR data are not used."))
}
if (length(cell_ids) == 0) {
    stop(
        "No TCR data found in the Seurat object. ",
        "Please use scRepertiore::combineExpression() to generate the Seurat object with TCR data."
    )
}
expr <- as.matrix(expr)[, tcrdata$contig_id, drop=FALSE]

# Write input files
log$info("Writing input files ...")
write.table(tcrdata, file.path(tessa_dir, "tcrdata.txt"), sep=",", quote=FALSE, row.names=FALSE)
write.table(expr, file.path(tessa_dir, "exprdata.txt"), sep=",", quote=FALSE, row.names=TRUE, col.names=TRUE)

### End preparing input files for TESSA

### Start running TESSA
log$info("Running TESSA ...")

# The original TESSA uses a python wrapper to run the encoder and tessa model
# here we run those two steps directly here

log$info("- Running encoder ...")
cmd_encoder <- paste(
    python,
    file.path(tessa_srcdir, "BriseisEncoder.py"),
    "-tcr",
    file.path(tessa_dir, "tcrdata.txt"),
    "-model",
    file.path(tessa_srcdir, "TrainedEncoder.h5"),
    "-embeding_vectors",
    file.path(tessa_srcdir, "Atchley_factors.csv"),
    "-output_TCR",
    file.path(tessa_dir, "tcr_encoded.txt"),
    "-output_log",
    file.path(tessa_dir, "tcr_encoder.log")
)
cmd_encoder <- paste(
    cmd_encoder,
    "-output_VJ",
    file.path(tessa_dir, "tcr_vj.txt")
)

print("Running:")
print(cmd_encoder)
log$debug(paste("- ", cmd_encoder))

rc <- system(cmd_encoder)
if (rc != 0) {
    stop("Error: Failed to run encoder.")
}

log$info("- Running TESSA model ...")
source(file.path(tessa_srcdir, "real_data.R"))

tessa <- run_tessa(
    tessa_srcdir,
    file.path(tessa_dir, "exprdata.txt"),
    file.path(tessa_dir, "tcr_encoded.txt"),
    file.path(tessa_dir, "tcrdata.txt"),
    result_dir,
    within_sample,
    (if (!predefined_b) NULL else file.path(tessa_srcdir, "fixed_b.csv")),
    max_iter = max_iter
)

# Save TESSA results
log$info("Saving TESSA results ...")
cells <- rownames(sobj@meta.data)
sobj@meta.data <- sobj@meta.data %>%
    mutate(
        TESSA_Cluster = tessa$meta[
            match(cells, tessa$meta$barcode),
            "cluster_number"
        ]
    ) %>%
    add_count(TESSA_Cluster, name = "TESSA_Cluster_Size")
rownames(sobj@meta.data) <- cells

if (save_tessa) {
    sobj@misc$tessa <- tessa
}
save_obj(sobj, outfile)

# Post analysis
log$info("Post analysis ...")
plot_tessa(tessa, result_dir)
plot_Tessa_clusters(tessa, result_dir)

p <- tessa$meta %>%
    dplyr::select(barcode, TESSA_Cluster = cluster_number) %>%
    add_count(TESSA_Cluster, name = "TESSA_Cluster_Size") %>%
    plotthis::Histogram(x = "TESSA_Cluster_Size")

res <- 100
height <- attr(p, "height") * res
width <- attr(p, "width") * res
prefix <- file.path(result_dir, "Cluster_size_dist")
save_plot(p, prefix, devpars = list(width = width, height = height, res = res))

reporter$add(
    list(
        src = file.path(result_dir, "Cluster_size_dist.png"),
        descr = "Histogram of cluster size distribution",
        download = file.path(result_dir, "Cluster_size_dist.pdf")
    ),
    list(
        src = file.path(result_dir, "clone_size.png"),
        descr = "Center cluster size vs. non-center cluster size"
    ),
    list(
        src = file.path(result_dir, "exp_TCR_pair_plot.png"),
        descr = "Expression-TCR distance plot"
    ),
    list(
        src = file.path(result_dir, "TCR_dist_density.png"),
        descr = "TCR distance density plot"
    ),
    list(
        src = file.path(result_dir, "TCR_explore.png"),
        descr = "Exploratory plot at the TCR level"
    ),
    list(
        src = file.path(result_dir, "TCR_explore_clusters.png"),
        descr = "TESSA clusters"
    ),
    h1 = "TESSA Results",
    ui = "table_of_images"
)

reporter$save(joboutdir)
