{% from "utils/gsea.liq" import fgsea_report -%}
{% from "utils/misc.liq" import report_jobs -%}
<script>
    import { Image, DataTable } from "$libs";
</script>

{%- macro report_job(job, h=1) -%}
{% for secdir in job.out.outdir | glob: "*" %}
    {%- set section = secdir | basename -%}
    {%- if section != "DEFAULT" -%}
    <h{{h}}>{{section}}</h{{h}}>
    {%- else -%}
    {%- set h = h - 1 -%}
    {%- endif -%}
    {% for casedir in secdir | glob: "*" %}
        {%- set case = casedir | basename -%}
        <h{{h+1}}>{{case}}</h{{h+1}}>
        {{ fgsea_report(casedir, h + 2) }}
    {%- endfor -%}
{% endfor %}
{%- endmacro -%}

{%- macro head_job(job) -%}
<h1>{{job.in.srtobj | stem0 | escape}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
