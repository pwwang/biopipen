source("{{biopipen_dir}}/utils/misc.R")

ngenes <- {{in.ngenes | r}}
nsamples <- {{in.nsamples | r}}
outdir <- {{out.outdir | r}}
seed <- {{envs.seed | r}}
ncores <- {{envs.ncores | r}}

{% if envs.tool.lower()  == "ruvcorr" %}
{% include biopipen_dir + "/scripts/rnaseq/Simulation-RUVcorr.R" %}
{% elif envs.tool.lower() == "esco" %}
{% include biopipen_dir + "/scripts/rnaseq/Simulation-ESCO.R" %}
{% else %}
stop("Unknown tool: {{envs.tool}}, only 'RUVcorr' and 'ESCO' are supported.")
{% endif %}
