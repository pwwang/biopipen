{% if envs.tool == "rmagic" %}
{% include biopipen_dir + "/scripts/scrna/ExprImpute-rmagic.R" %}
{% elif envs.tool == "scimpute" %}
{% include biopipen_dir + "/scripts/scrna/ExprImpute-scimpute.R" %}
{% endif %}
