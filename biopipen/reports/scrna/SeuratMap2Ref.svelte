{% from "utils/misc.liq" import report_jobs, table_of_images -%}
<script>
    import { Image, DataTable } from "$libs";
</script>

{%- macro report_job(job, h=1) -%}

<h{{h}}>UMAPs</h{{h}}>
{% set imgs = [] %}
{% set caps = [] %}
{% for png in job.outdir | glob: "UMAPs-*.png" %}
    {% set pdf = png | regex_replace: "\\.png$", ".pdf" %}
    {% set stm = png | stem %}
    {% set _ = imgs.append({"src": png, "download": pdf}) %}
    {% set _ = caps.append(stm | replace: "UMAPs-", "") %}
{% endfor %}
{{ table_of_images(imgs, caps) }}

<h{{h}}>Mapping Score</h{{h}}>
<Image
    src="{{job.outdir | joinpath: 'mapping_score.png'}}"
    download="{{job.outdir | joinpath: 'mapping_score.pdf'}}"
    />

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
