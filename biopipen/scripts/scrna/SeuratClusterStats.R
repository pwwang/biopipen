library(Seurat)
library(rlang)
library(dplyr)
library(tidyr)
library(tibble)
library(glue)
library(forcats)
library(tidyseurat)
library(gglogger)
library(scplotter)
library(biopipen.utils)

log <- get_logger()
reporter <- get_reporter()

srtfile <- {{in.srtobj | r}}
outdir <- {{out.outdir | r}}
joboutdir <- {{job.outdir | r}}
mutaters <- {{envs.mutaters | r}}
cache <- {{envs.cache | r}}

if (isTRUE(cache)) { cache = joboutdir }

log$info("Loading Seurat object ...")
srtobj = read_obj(srtfile)

log$info("Applying mutaters if any ...")
if (!is.null(mutaters) && length(mutaters) > 0) {
    srtobj@meta.data = srtobj@meta.data %>%
        mutate(!!!lapply(mutaters, parse_expr))
}

############## clustree ##############
clustrees_defaults <- {{envs.clustrees_defaults | r: todot="-"}}
clustrees <- {{envs.clustrees | r: todot="-", skip=1}}
{% include biopipen_dir + "/scripts/scrna/SeuratClusterStats-clustree.R" %}

############## stats ##############
stats_defaults = {{envs.stats_defaults | r: todot="-"}}
stats = {{envs.stats | r: todot="-", skip=1}}
{% include biopipen_dir + "/scripts/scrna/SeuratClusterStats-stats.R" %}

############## ngenes ##############
ngenes_defaults <- {{envs.ngenes_defaults | r: todot="-"}}
ngenes <- {{envs.ngenes | r: todot="-", skip=1}}
{% include biopipen_dir + "/scripts/scrna/SeuratClusterStats-ngenes.R" %}

############## features ##############
features_defaults = {{envs.features_defaults | r: todot="-"}}
features = {{envs.features | r: todot="-", skip=1}}
{% include biopipen_dir + "/scripts/scrna/SeuratClusterStats-features.R" %}

############## dimplots ##############
dimplots_defaults = {{envs.dimplots_defaults | r: todot="-"}}
dimplots = {{envs.dimplots | r: todot="-", skip=1}}
{% include biopipen_dir + "/scripts/scrna/SeuratClusterStats-dimplots.R" %}

reporter$save(joboutdir)
