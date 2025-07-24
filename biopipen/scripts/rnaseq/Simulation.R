ngenes <- {{in.ngenes | r}}
nsamples <- {{in.nsamples | r}}
outfile <- {{out.outfile | r}}
outdir <- {{out.outdir | r}}
seed <- {{envs.seed | r}}
ncores <- {{envs.ncores | r}}
transpose_output <- {{envs.transpose_output | r}}
index_start <- {{envs.index_start | r}}

{% if envs.tool.lower()  == "ruvcorr" %}
{% include biopipen_dir + "/scripts/rnaseq/Simulation-RUVcorr.R" %}
{% elif envs.tool.lower() == "esco" %}
{% include biopipen_dir + "/scripts/rnaseq/Simulation-ESCO.R" %}
{% else %}
stop("Unknown tool: {{envs.tool}}, only 'RUVcorr' and 'ESCO' are supported.")
{% endif %}

colnames(simulated) <- paste0("Sample", index_start + 0:(nsamples - 1))
if (transpose_output) { simulated <- t(simulated) }

write.table(simulated, file = outfile, sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)
