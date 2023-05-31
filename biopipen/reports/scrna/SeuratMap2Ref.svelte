{% from "utils/misc.liq" import report_jobs -%}
<script>
    import { Image } from "$libs";
</script>

{%- macro report_job(job, h=1) -%}
<h{{h}}>Reference UMAP</h{{h}}>
{% for refumap in job.outdir | glob: "Reference_UMAP_*.png" %}
    <Image src="{{refumap}}" />
{% endfor %}

<h{{h}}>Query UMAP</h{{h}}>
{% for qryumap in job.outdir | glob: "Query_UMAP_*.png" %}
    <Image src="{{qryumap}}" />
{% endfor %}
{%- endmacro -%}

{%- macro head_job(job) -%}
<h1>{{job.in.sobjfile | stem0 | escape}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
