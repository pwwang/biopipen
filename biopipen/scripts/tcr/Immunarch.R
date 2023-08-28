source("{{biopipen_dir}}/utils/misc.R")
# Basic analysis and clonality
# TODO: How about TRA chain?
library(rlang)
library(immunarch)  # 0.8.0 or 0.9.0 tested
library(ggplot2)
library(ggprism)
library(dplyr)
library(tidyr)
library(tibble)

theme_set(theme_prism())

immfile = {{ in.immdata | quote }}
outdir = {{ out.outdir | quote }}
mutaters = {{ envs.mutaters | r }}

immdata = readRDS(immfile)
if (!is.null(mutaters) && length(mutaters) > 0) {
    immdata$meta = mutate(immdata$meta, !!!lapply(mutaters, parse_expr))
}

n_samples = length(immdata$data)

# Basic analysis
{% include biopipen_dir + "/scripts/tcr/Immunarch-basic.R" %}

# Clonality
{% include biopipen_dir + "/scripts/tcr/Immunarch-clonality.R" %}

# Overlap
{% include biopipen_dir + "/scripts/tcr/Immunarch-overlap.R" %}

# Gene usage
{% include biopipen_dir + "/scripts/tcr/Immunarch-geneusage.R" %}

# Spectratyping
{% include biopipen_dir + "/scripts/tcr/Immunarch-spectratyping.R" %}

# Diversity estimation
{% include biopipen_dir + "/scripts/tcr/Immunarch-diversity.R" %}

# Clonotype tracking
{% include biopipen_dir + "/scripts/tcr/Immunarch-tracking.R" %}

# K-mer analysis
{% include biopipen_dir + "/scripts/tcr/Immunarch-kmer.R" %}
