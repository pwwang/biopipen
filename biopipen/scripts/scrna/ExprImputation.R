{% if envs.tool == "rmagic" %}
{% include biopipen_dir + "/scripts/scrna/ExprImputation-rmagic.R" %}
{% elif envs.tool == "scimpute" %}
{% include biopipen_dir + "/scripts/scrna/ExprImputation-scimpute.R" %}
{% elif envs.tool == "alra" %}
{% include biopipen_dir + "/scripts/scrna/ExprImputation-alra.R" %}
{% endif %}
