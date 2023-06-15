{% from "utils/misc.liq" import report_jobs, table_of_images -%}
<script>
    import { Image } from "$libs";
</script>

{%- macro report_job(job, h=1) -%}
<h{{h}}>Reference UMAP</h{{h}}>
{% set imgs = job.outdir | glob: "Reference_UMAP_*.png" %}
{{ table_of_images(imgs) }}

<h{{h}}>Query UMAP</h{{h}}>
{% set imgs = job.outdir | glob: "Query_UMAP_*.png" %}
{{ table_of_images(imgs) }}
{%- endmacro -%}

{%- macro head_job(job) -%}
<h1>{{job.in.sobjfile | stem0 | escape}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
