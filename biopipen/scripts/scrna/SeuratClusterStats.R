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

############## clustree ##############
clustrees_defaults <- {{envs.clustrees_defaults | r}}
clustrees <- {{envs.clustrees | r}}
{% set sourcefile = biopipen_dir | joinpaths: "scripts", "scrna", "SeuratClusterStats-clustree.R" %}
# force update when sourcefile changes
# {{ sourcefile | getmtime }}
source("{{sourcefile}}")

############## stats ##############
stats_defaults = {{envs.stats_defaults | r: todot="-"}}
stats = {{envs.stats | r: todot="-", skip=1}}
{% set sourcefile = biopipen_dir | joinpaths: "scripts", "scrna", "SeuratClusterStats-stats.R" %}
# {{ sourcefile | getmtime }}
source("{{sourcefile}}")

############## hists ##############
hists_defaults <- {{envs.hists_defaults | r: todot="-"}}
hists <- {{envs.hists | r: todot="-", skip=1}}
{% set sourcefile = biopipen_dir | joinpaths: "scripts", "scrna", "SeuratClusterStats-hists.R" %}
# {{ sourcefile | getmtime }}
source("{{sourcefile}}")

############## ngenes ##############
ngenes_defaults <- {{envs.ngenes_defaults | r: todot="-"}}
ngenes <- {{envs.ngenes | r: todot="-", skip=1}}
{% set sourcefile = biopipen_dir | joinpaths: "scripts", "scrna", "SeuratClusterStats-ngenes.R" %}
# {{ sourcefile | getmtime }}
source("{{sourcefile}}")

############## features ##############
features_defaults = {{envs.features_defaults | r: todot="-"}}
features = {{envs.features | r: todot="-", skip=1}}
{% set sourcefile = biopipen_dir | joinpaths: "scripts", "scrna", "SeuratClusterStats-features.R" %}
# {{ sourcefile | getmtime }}
source("{{sourcefile}}")

############## dimplots ##############
dimplots_defaults = {{envs.dimplots_defaults | r: todot="-"}}
dimplots = {{envs.dimplots | r: todot="-", skip=1}}
{% set sourcefile = biopipen_dir | joinpaths: "scripts", "scrna", "SeuratClusterStats-dimplots.R" %}
# {{ sourcefile | getmtime }}
source("{{sourcefile}}")

save_report(joboutdir)
