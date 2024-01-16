source("{{biopipen_dir}}/utils/misc.R")
source("{{biopipen_dir}}/utils/mutate_helpers.R")
source("{{biopipen_dir}}/utils/plot.R")
library(Seurat)
library(rlang)
library(dplyr)
library(tidyr)
library(tibble)
library(ggprism)
library(ggrepel)
library(tidyseurat)
library(circlize)
library(ComplexHeatmap)

srtfile <- {{in.srtobj | r}}
outdir <- {{out.outdir | r}}
joboutdir <- {{job.outdir | r}}
mutaters <- {{envs.mutaters | r}}

log_info("Loading Seurat object ...")
srtobj = readRDS(srtfile)

log_info("Applying mutaters if any ...")
if (!is.null(mutaters) && length(mutaters) > 0) {
    srtobj@meta.data = srtobj@meta.data %>%
        mutate(!!!lapply(mutaters, parse_expr))
}

############## stats ##############
{% include biopipen_dir + "/scripts/scrna/SeuratClusterStats-stats.R" %}

############## hists ##############
{% include biopipen_dir + "/scripts/scrna/SeuratClusterStats-hists.R" %}

############## ngenes ##############
{% include biopipen_dir + "/scripts/scrna/SeuratClusterStats-ngenes.R" %}

############## features ##############
{% include biopipen_dir + "/scripts/scrna/SeuratClusterStats-features.R" %}

############## dimplots ##############
{% include biopipen_dir + "/scripts/scrna/SeuratClusterStats-dimplots.R" %}

save_report(joboutdir)
