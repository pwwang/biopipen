{% if envs.tool == "rmagic" %}
{% include biopipen_dir + "/scripts/scrna/ExprImpution-rmagic.R" %}
{% elif envs.tool == "scimpute" %}
{% include biopipen_dir + "/scripts/scrna/ExprImpution-scimpute.R" %}
{% elif envs.tool == "alra" %}
{% include biopipen_dir + "/scripts/scrna/ExprImpution-alra.R" %}
{% endif %}
