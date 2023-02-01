{% if envs.tool == "rmagic" %}
{% include biopipen_dir + "/scripts/scrna/ExprImpute-rmagic.R" %}
{% elif envs.tool == "scimpute" %}
{% include biopipen_dir + "/scripts/scrna/ExprImpute-scimpute.R" %}
{% elif envs.tool == "alra" %}
{% include biopipen_dir + "/scripts/scrna/ExprImpute-alra.R" %}
{% endif %}
