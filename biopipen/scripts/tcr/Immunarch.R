source("{{biopipen_dir}}/utils/misc.R")
source("{{biopipen_dir}}/utils/single_cell.R")
# Basic analysis and clonality
# TODO: How about TRA chain?
library(rlang)
library(immunarch)  # 0.8.0 or 0.9.0 tested
library(ggplot2)
library(patchwork)
library(ggprism)
library(dplyr)
library(glue)
library(tidyr)
library(tibble)
library(logger)

log_info("Loading arguments ...")
theme_set(theme_prism())

immfile = {{ in.immdata | r }}
metafile = {{ in.metafile | r }}
outdir = {{ out.outdir | r }}
joboutdir = {{ job.outdir | r }}
mutaters = {{ envs.mutaters | r }}
prefix = {{ envs.prefix | r }}

log_info("Loading immdata ...")
immdata = readRDS(immfile)

if (is.null(prefix)) { prefix = immdata$prefix }
if (is.null(prefix)) { prefix = "" }

log_info("Expanding immdata ...")
exdata = expand_immdata(immdata)

log_info("Loading metadata if provided ...")
if (!is.null(metafile) && (endsWith(metafile, ".rds") || endsWith(metafile, ".RDS"))) {
    meta = readRDS(metafile)@meta.data
} else if (!is.null(metafile) && nchar(metafile) > 0) {
    meta = read.table(metafile, sep = "\t", header = TRUE, row.names = 1)
} else {
    meta = NULL
}

log_info("Merging metadata if provided ...")
if (!is.null(meta)) {
    cell_ids = glue_data(exdata, paste0(prefix, "{Barcode}"))
    dup_names = intersect(colnames(meta), colnames(exdata))
    if (length(dup_names) > 0) {
        warning(paste0("Duplicated column names in meta: ", paste(dup_names, collapse = ", ")), immediate. = TRUE)
        meta = rename(meta, !!!setNames(dup_names, paste0(dup_names, "_meta")))
    }
    exdata = cbind(exdata, meta[cell_ids, , drop = FALSE])
    rm(cell_ids)
}
rm(meta)

log_info("Mutating data if `envs.mutaters` is provided ...")
if (!is.null(mutaters) && length(mutaters) > 0) {
    exdata = mutate(exdata, !!!lapply(mutaters, parse_expr))
    immdata = immdata_from_expanded(exdata)
}

n_samples = length(immdata$data)

##################
# Basic analysis #
##################
{% include biopipen_dir + "/scripts/tcr/Immunarch-basic.R" %}

##################
# Clonality      #
##################
{% include biopipen_dir + "/scripts/tcr/Immunarch-clonality.R" %}

##################
# Overlap        #
##################
{% include biopipen_dir + "/scripts/tcr/Immunarch-overlap.R" %}

##################
# Gene usage     #
##################
{% include biopipen_dir + "/scripts/tcr/Immunarch-geneusage.R" %}

##################
# Spectratyping  #
##################
{% include biopipen_dir + "/scripts/tcr/Immunarch-spectratyping.R" %}

########################
# Diversity estimation #
########################
{% include biopipen_dir + "/scripts/tcr/Immunarch-diversity.R" %}

######################
# Clonotype tracking #
######################
{% include biopipen_dir + "/scripts/tcr/Immunarch-tracking.R" %}

######################
# K-mer analysis     #
######################
{% include biopipen_dir + "/scripts/tcr/Immunarch-kmer.R" %}

######################
# VJ junction        #
######################
{% include biopipen_dir + "/scripts/tcr/Immunarch-vjjunc.R" %}

save_report(joboutdir)
