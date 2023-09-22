set.seed(8525)

{% if envs.tool == "hitype" %}
{% include biopipen_dir + "/scripts/scrna/CellTypeAnnotation-hitype.R" %}
{% elif envs.tool == "sctype" %}
{% include biopipen_dir + "/scripts/scrna/CellTypeAnnotation-sctype.R" %}
{% elif envs.tool == "sccatch" %}
{% include biopipen_dir + "/scripts/scrna/CellTypeAnnotation-sccatch.R" %}
{% elif envs.tool == "direct" %}
{% include biopipen_dir + "/scripts/scrna/CellTypeAnnotation-direct.R" %}
{% else %}
stop(paste0("Unknown tool: ", {{envs.tool}}))
{% endif %}
