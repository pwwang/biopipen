source("{{biopipen_dir}}/utils/misc.R")
source("{{biopipen_dir}}/utils/plot.R")
library(slugify)
library(Seurat)
library(rlang)
library(dplyr)
library(tibble)
library(ggprism)
library(ggrepel)
library(tidyseurat)

srtfile = {{in.srtobj | r}}
outdir = {{out.outdir | r}}
joboutdir = {{job.outdir | r}}

log_info("Loading Seurat object ...")
srtobj = readRDS(srtfile)

{% include biopipen_dir + "/scripts/scrna/SeuratClusterStats-stats.R" %}
{% include biopipen_dir + "/scripts/scrna/SeuratClusterStats-features.R" %}
{% include biopipen_dir + "/scripts/scrna/SeuratClusterStats-dimplots.R" %}

save_report(joboutdir)
