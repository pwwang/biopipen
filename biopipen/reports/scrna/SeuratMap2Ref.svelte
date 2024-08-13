{% from "utils/misc.liq" import report_jobs, table_of_images -%}
<script>
    import { Image, DataTable } from "$libs";
</script>

{%- macro report_job(job, h=1) -%}

<h{{h}}>UMAPs</h{{h}}>
{% set imgs = job.outdir | glob: "UMAPs-*.png" %}
{{ table_of_images(imgs) }}

<h{{h}}>Stats</h{{h}}>
{% for stfile in job.outdir | glob: "stats-*.txt" %}
    <h{{h+1}}>{{stfile | stem | replace: "stats-", ""}}</h{{h+1}}>
    <DataTable src="{{stfile}}" data={ {{stfile | datatable: sep="\t"}} } />
{% endfor %}

{%- endmacro -%}

{%- macro head_job(job) -%}
<h1>{{job.in.sobjfile | stem0 | escape}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
