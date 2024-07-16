{{ biopipen_dir | joinpaths: "utils", "misc.R" | source_r }}
{{ biopipen_dir | joinpaths: "utils", "mutate_helpers.R" | source_r }}
{{ biopipen_dir | joinpaths: "utils", "plot.R" | source_r }}

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
{{ biopipen_dir | joinpaths: "scripts", "scrna", "SeuratClusterStats-clustree.R" | source_r }}

############## stats ##############
stats_defaults = {{envs.stats_defaults | r: todot="-"}}
stats = {{envs.stats | r: todot="-", skip=1}}
{{ biopipen_dir | joinpaths: "scripts", "scrna", "SeuratClusterStats-stats.R" | source_r }}

############## hists ##############
hists_defaults <- {{envs.hists_defaults | r: todot="-"}}
hists <- {{envs.hists | r: todot="-", skip=1}}
{{ biopipen_dir | joinpaths: "scripts", "scrna", "SeuratClusterStats-hists.R" | source_r }}

############## ngenes ##############
ngenes_defaults <- {{envs.ngenes_defaults | r: todot="-"}}
ngenes <- {{envs.ngenes | r: todot="-", skip=1}}
{{ biopipen_dir | joinpaths: "scripts", "scrna", "SeuratClusterStats-ngenes.R" | source_r }}

############## features ##############
features_defaults = {{envs.features_defaults | r: todot="-"}}
features = {{envs.features | r: todot="-", skip=1}}
{{ biopipen_dir | joinpaths: "scripts", "scrna", "SeuratClusterStats-features.R" | source_r }}

############## dimplots ##############
dimplots_defaults = {{envs.dimplots_defaults | r: todot="-"}}
dimplots = {{envs.dimplots | r: todot="-", skip=1}}
{{ biopipen_dir | joinpaths: "scripts", "scrna", "SeuratClusterStats-dimplots.R" | source_r }}

save_report(joboutdir)
