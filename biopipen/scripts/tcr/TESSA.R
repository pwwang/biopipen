source("{{biopipen_dir}}/utils/misc.R")
source("{{biopipen_dir}}/utils/single_cell.R")

library(glue)
library(dplyr)
library(tidyr)
library(tibble)
library(immunarch)
library(Seurat)
library(ggplot2)
library(ggprism)

immfile <- {{in.immdata | r}}
exprfile <- {{in.srtobj | r}}
outfile <- {{out.outfile | r}}
joboutdir <- {{job.outdir | r}}
python <- {{envs.python | r}}
prefix <- {{envs.prefix | r}}
within_sample <- {{envs.within_sample | r}}
assay <- {{envs.assay | r}}
predefined_b <- {{envs.predefined_b | r}}
max_iter <- {{envs.max_iter | int}}
save_tessa <- {{envs.save_tessa | r}}
tessa_srcdir <- "{{biopipen_dir}}/scripts/tcr/TESSA_source"

outdir <- dirname(outfile)
result_dir <- file.path(outdir, "result")
tessa_dir <- file.path(outdir, "tessa")
if (!dir.exists(result_dir)) dir.create(result_dir)
if (!dir.exists(tessa_dir)) dir.create(tessa_dir)

### Start preparing input files for TESSA
# Prepare input files
log_info("Preparing TCR input file ...")
# If immfile endswith .rds, then it is an immunarch object
if (endsWith(tolower(immfile), ".rds")) {
    immdata <- readRDS(immfile)
    if (is.null(prefix)) { prefix = immdata$prefix }
    if (is.null(prefix)) { prefix = "" }
    tcrdata <- expand_immdata(immdata) %>%
        mutate(Barcode = glue(paste0(prefix, "{Barcode}")))
    rm(immdata)
} else {
    tcrdata <- read.table(immfile, sep="\t", header=TRUE, row.names=1) %>%
        rownames_to_column("Barcode")
}

has_VJ <- "V.name" %in% colnames(tcrdata) && "J.name" %in% colnames(tcrdata)

if (has_VJ) {
    tcrdata <- tcrdata %>% dplyr::mutate(
        v_gene = sub("-\\d+$", "", V.name),
        j_gene = sub("-\\d+$", "", J.name)
    ) %>% dplyr::select(
        contig_id = Barcode,
        cdr3 = CDR3.aa,
        v_gene,
        j_gene,
        sample = Sample
    )
} else {
    tcrdata <- tcrdata %>% dplyr::select(
        contig_id = Barcode,
        cdr3 = CDR3.aa,
        sample = Sample
    )
}


log_info("Preparing expression input file ...")
is_seurat <- endsWith(tolower(exprfile), ".rds")
is_gz <- endsWith(tolower(exprfile), ".gz")

if (is_seurat) {
    sobj <- readRDS(exprfile)
    expr <- GetAssayData(sobj, layer = "data")
} else if (is_gz) {
    expr <- read.table(gzfile(exprfile), sep="\t", header=TRUE, row.names=1)
} else {
    expr <- read.table(exprfile, sep="\t", header=TRUE, row.names=1)
}

cell_ids <- intersect(tcrdata$contig_id, colnames(expr))
# Warning about unused cells
unused_tcr_cells <- setdiff(tcrdata$contig_id, cell_ids)
unused_expr_cells <- setdiff(colnames(expr), cell_ids)
if (length(unused_tcr_cells) > 0) {
    log_warn(glue("{length(unused_tcr_cells)}/{nrow(tcrdata)} TCR cells are not used."))
}
if (length(unused_expr_cells) > 0) {
    log_warn(glue("{length(unused_expr_cells)}/{ncol(expr)} expression cells are not used."))
}
if (length(cell_ids) == 0) {
    stop(paste0(
        "No common cells between TCR and expression data. ",
        "Are you using the correct `envs.prefix` here or in `ImmunarchLoading`?"
    ))
}
tcrdata <- tcrdata[tcrdata$contig_id %in% cell_ids, , drop=FALSE]
expr <- as.matrix(expr)[, tcrdata$contig_id, drop=FALSE]

# Write input files
log_info("Writing input files ...")
write.table(tcrdata, file.path(tessa_dir, "tcrdata.txt"), sep=",", quote=FALSE, row.names=FALSE)
write.table(expr, file.path(tessa_dir, "exprdata.txt"), sep=",", quote=FALSE, row.names=TRUE, col.names=TRUE)

### End preparing input files for TESSA

### Start running TESSA
log_info("Running TESSA ...")

# The original TESSA uses a python wrapper to run the encoder and tessa model
# here we run those two steps directly here

log_info("- Running encoder ...")
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
if (has_VJ) {
    cmd_encoder <- paste(
        cmd_encoder,
        "-output_VJ",
        file.path(tessa_dir, "tcr_vj.txt")
    )
}
print("Running:")
print(cmd_encoder)
log_debug(paste("- ", cmd_encoder))

rc <- system(cmd_encoder)
if (rc != 0) {
    stop("Error: Failed to run encoder.")
}

log_info("- Running TESSA model ...")
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
log_info("Saving TESSA results ...")
if (is_seurat) {
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
    saveRDS(sobj, outfile)
} else {
    out <- tessa$meta %>%
        dplyr::select(barcode, TESSA_Cluster = cluster_number) %>%
        add_count(TESSA_Cluster, name = "TESSA_Cluster_Size")
    write.table(out, outfile, sep="\t", quote=FALSE, row.names=FALSE, col.names=TRUE)
}

# Post analysis
log_info("Post analysis ...")
plot_tessa(tessa, result_dir)
plot_Tessa_clusters(tessa, result_dir)

p <- tessa$meta %>%
    dplyr::select(barcode, TESSA_Cluster = cluster_number) %>%
    add_count(TESSA_Cluster, name = "TESSA_Cluster_Size") %>%
    ggplot(aes(x = TESSA_Cluster_Size)) +
    geom_histogram(binwidth = 1) +
    theme_prism()

png(file.path(result_dir, "Cluster_size_dist.png"), width=8, height=8, units="in", res=100)
print(p)
dev.off()

add_report(
    list(
        src = file.path(result_dir, "Cluster_size_dist.png"),
        descr = "Histogram of cluster size distribution"
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

save_report(joboutdir)
