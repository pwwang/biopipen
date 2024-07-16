{{ biopipen_dir | joinpaths: "utils", "misc.R" | source_r }}
{{ biopipen_dir | joinpaths: "utils", "single_cell.R" | source_r }}
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
volumes = {{envs.volumes | r}}
lens = {{envs.lens | r}}
counts = {{envs.counts | r}}
{{ biopipen_dir | joinpaths: "scripts", "tcr", "Immunarch-basic.R" | source_r }}

##################
# Clonality      #
##################
top_clones = {{envs.top_clones | r}}
rare_clones = {{envs.rare_clones | r}}
hom_clones = {{envs.hom_clones | r}}
{{ biopipen_dir | joinpaths: "scripts", "tcr", "Immunarch-clonality.R" | source_r }}

##################
# Overlap        #
##################
overlaps = {{ envs.overlaps | r: todot="-" }}
{{ biopipen_dir | joinpaths: "scripts", "tcr", "Immunarch-overlap.R" | source_r }}

##################
# Gene usage     #
##################
gene_usages = {{ envs.gene_usages | r: todot="-" }}
{{ biopipen_dir | joinpaths: "scripts", "tcr", "Immunarch-geneusage.R" | source_r }}

##################
# Spectratyping  #
##################
spects = {{ envs.spects | r }}
{{ biopipen_dir | joinpaths: "scripts", "tcr", "Immunarch-spectratyping.R" | source_r }}

########################
# Diversity estimation #
########################
div_method = {{envs.divs.method | default: "gini" | r}}
div_by = {{envs.divs.by | default: None | r}}
div_plot_type = {{envs.divs.plot_type | default: "bar" | r}}
div_order = {{envs.divs.order | default: [] | r}}
div_args = {{envs.divs.args | default: {} | r: todot="-"}}
div_test = {{envs.divs.test | default: None | r}}
div_cases = {{envs.divs.cases | default: {} | r: todot="-"}}
div_devpars = {{envs.divs.devpars | default: None | r}}
div_separate_by = {{envs.divs.separate_by | default: None | r}}
div_split_by = {{envs.divs.split_by | default: None | r}}
div_split_order = {{envs.divs.split_order | default: None | r}}
div_align_x = {{envs.divs.align_x | default: False | r}}
div_align_y = {{envs.divs.align_y | default: False | r}}
div_subset = {{envs.divs.subset | default: None | r}}
div_log = {{envs.divs.log | default: False | r}}
div_ncol = {{envs.divs.ncol | default: 2 | r}}
div_ymin = {{envs.divs.ymin | default: None | r}}
div_ymax = {{envs.divs.ymax | default: None | r}}

{{ biopipen_dir | joinpaths: "scripts", "tcr", "immunarch-patched.R" | source_r }}
{{ biopipen_dir | joinpaths: "scripts", "tcr", "Immunarch-diversity.R" | source_r }}

######################
# Clonotype tracking #
######################
trackings = {{ envs.trackings | r }}
{{ biopipen_dir | joinpaths: "scripts", "tcr", "Immunarch-tracking.R" | source_r }}

######################
# K-mer analysis     #
######################
kmers = {{ envs.kmers | r: todot="-" }}
{{ biopipen_dir | joinpaths: "scripts", "tcr", "Immunarch-kmer.R" | source_r }}

######################
# VJ junction        #
######################
vj_juncs <- {{envs.vj_junc | r}}
{{ biopipen_dir | joinpaths: "scripts", "tcr", "Immunarch-vjjunc.R" | source_r }}

save_report(joboutdir)
