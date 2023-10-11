source("{{biopipen_dir}}/utils/misc.R")

library(glue)
library(dplyr)
library(tidyr)
library(immunarch)
library(Seurat)
library(ggplot2)
library(ggprism)

immfile <- {{in.immdata | r}}
exprfile <- {{in.srtobj | r}}
outfile <- {{out.outfile | r}}
python <- {{envs.python | r}}
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
print("Preparing TCR input file ...")
immdata <- readRDS(immfile)

has_VJ <- "V.name" %in% colnames(immdata$data[[1]]) && "J.name" %in% colnames(immdata$data[[1]])
# Merge all samples
tcrdata <- do_call(rbind, lapply(seq_len(nrow(immdata$meta)), function(i) {
    # Clones  Proportion   CDR3.aa                       Barcode
    # 5      4 0.008583691 CAVRDTGNTPLVF;CASSEYSNQPQHF   GTTCGGGCACTTACGA-1;TCTCTAAGTACCAGTT-1
    # 6      4 0.008583691 CALTQAAGNKLTF;CASRPEDLRGQPQHF GCTTGAAGTCGGCACT-1;TACTCGCTCCTAAGTG-1
    if (has_VJ) {
        cldata = immdata$data[[i]][, c("Barcode", "CDR3.aa", "V.name", "J.name")]
    } else {
        cldata = immdata$data[[i]][, c("Barcode", "CDR3.aa")]
    }
    # # A tibble: 4 × 5
    # Sample                  Patient     Timepoint Tissue
    # <chr>                   <chr>       <chr>     <chr>
    # 1 MC1685Pt011-Baseline-PB MC1685Pt011 Baseline  PB
    mdata = as.list(immdata$meta[i, , drop=FALSE])
    for (mname in names(mdata)) {
        assign(mname, mdata[[mname]])
    }

    cldata %>%
        separate_rows(Barcode, sep=";") %>%
        # Just in case there are duplicated barcodes
        distinct(Barcode, .keep_all = TRUE) %>%
        mutate(Barcode = glue("{{envs.prefix}}{Barcode}"), sample = Sample)
}))
if (has_VJ) {
    tcrdata <- tcrdata %>% dplyr::mutate(
        v_gene = sub("-\\d+$", "", V.name),
        j_gene = sub("-\\d+$", "", J.name)
    ) %>% dplyr::select(
        contig_id = Barcode,
        cdr3 = CDR3.aa,
        v_gene,
        j_gene,
        sample
    )
} else {
    tcrdata <- tcrdata %>% dplyr::select(
        contig_id = Barcode,
        cdr3 = CDR3.aa,
        sample
    )
}


print("Preparing expression input file ...")
is_seurat <- endsWith(tolower(exprfile), ".rds")
is_gz <- endsWith(tolower(exprfile), ".gz")

if (is_seurat) {
    sobj <- readRDS(exprfile)
    expr <- GetAssayData(sobj, slot = "data", assay = assay)
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
    warning(glue("{length(unused_tcr_cells)}/{nrow(tcrdata)} TCR cells are not used."), immediate. = TRUE)
}
if (length(unused_expr_cells) > 0) {
    warning(glue("{length(unused_expr_cells)}/{ncol(expr)} expression cells are not used."), immediate. = TRUE)
}
if (length(cell_ids) == 0) {
    stop("No common cells between TCR and expression data. Are you using the correct prefix?")
}
tcrdata <- tcrdata[tcrdata$contig_id %in% cell_ids, , drop=FALSE]
expr <- as.matrix(expr)[, tcrdata$contig_id, drop=FALSE]

# Write input files
print("Writing input files ...")
write.table(tcrdata, file.path(tessa_dir, "tcrdata.txt"), sep=",", quote=FALSE, row.names=FALSE)
write.table(expr, file.path(tessa_dir, "exprdata.txt"), sep=",", quote=FALSE, row.names=TRUE, col.names=TRUE)

### End preparing input files for TESSA

### Start running TESSA
print("Running TESSA ...")

# The original TESSA uses a python wrapper to run the encoder and tessa model
# here we run those two steps directly here

print("- Running encoder ...")
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
print(paste("- ", cmd_encoder))

rc <- system(cmd_encoder)
if (rc != 0) {
    stop("Error: Failed to run encoder.")
}

print("- Running TESSA model ...")
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
print("Saving TESSA results ...")
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
print("Post analysis ...")
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
