{% if envs.tool == "sctype" %}
{% include biopipen_dir + "/scripts/scrna/CellTypeAnnotate-sctype.R" %}
{% elif envs.tool == "direct" %}
{% include biopipen_dir + "/scripts/scrna/CellTypeAnnotate-direct.R" %}
{% else %}
stop(paste0("Unknown tool: ", {{envs.tool}}))
{% endif %}
