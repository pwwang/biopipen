{% from "utils/misc.liq" import table_of_images -%}
<script>
    import { Image } from "@@";
</script>

{% set images = [] %}
{% set names = [] %}

{% for job in jobs %}
{% set _ = images.append(job.out.outfile) %}
{% set config = job.in.configfile | read | toml_loads %}
{% set _ = names.append(config.name or "") %}
{% endfor %}

{{ table_of_images(images, names) }}
