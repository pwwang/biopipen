{% from "utils/misc.liq" import table_of_images, report_jobs -%}
<script>
    import { Image } from "$libs";
</script>

{%- macro report_job(job, h=1) -%}
{% set images = job.out.outdir | glob: "*.png" %}
{{ table_of_images(images) }}
{%- endmacro -%}

{%- macro head_job(job) -%}
{% if job.in.name %}
<h1>{{job.in.name | escape}}</h1>
{% else %}
<h1>Case {{job.index | plus: 1 | str}}</h1>
{% endif %}
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
