{% from "utils.liq" import report_jobs, table_of_images -%}
<script>
    import { Image } from "@@";
</script>

{%- macro report_job(job, h=1) -%}
{% for titlefile in job.out.outdir | joinpaths: "*.title" | glob %}
{%  assign idx = titlefile | stem %}
{%  assign png = job.out.outdir | joinpaths: idx + ".png" %}
<h{{h}}>{{ titlefile | read }}</h{{h}}>
<Image src={{png | quote}} />
{% endfor %}
{%- endmacro -%}

{%- macro head_job(job) -%}
<h1>{{job.out.outdir | stem | replace: ".exprs", ""}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
