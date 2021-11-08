{% from_ os import path %}
{% from "utils.liq" import report_jobs, table_of_images -%}

<script>
    import { Image } from "@@";
</script>

{%- macro report_job(job, h=1) -%}
{% for titlefile in job.out.outdir | joinpaths: "*.title" | glob %}
{%  assign idx = titlefile | stem %}
{%  assign boxplotpng = job.out.outdir | joinpaths: idx + "-boxplot.png" %}
{%  assign heatmappng = job.out.outdir | joinpaths: idx + "-heatmap.png" %}
<h{{h}}>{{ titlefile | read }}</h{{h}}>

{% if path.exists(boxplotpng) %}
<Image src={{boxplotpng | quote}} />
{% endif %}

{% if path.exists(heatmappng) %}
<Image src={{heatmappng | quote}} />
{% endif %}

{% endfor %}
{%- endmacro -%}

{%- macro head_job(job) -%}
<h1>{{job.out.outdir | stem | replace: ".exprs", ""}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
