{% from "utils/misc.liq" import report_jobs, table_of_images -%}
<script>
    import { Image } from "@@";
</script>

<h1>V-J usage plots</h1>

{% from_ os.path import basename %}
{% assign pngs = [] %}
{% assign caps = [] %}

{% for job in jobs %}
{% set _ = pngs.append(job.out.outfile) %}
{% set _ = caps.append(basename(job.out.outfile).replace(".fancyvj.wt.png", "")) %}
{% endfor %}

{{ table_of_images(pngs, caps) }}
