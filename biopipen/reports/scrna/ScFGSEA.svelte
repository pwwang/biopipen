{% from "utils/gsea.liq" import fgsea_report -%}
{% from "utils/misc.liq" import report_jobs -%}
<script>
    import { Image, DataTable } from "@@";
</script>

{%- macro report_job(job, h=1) -%}
{% for casedir in job.out.outdir | joinpaths: "*" | glob %}
{%  set case = casedir | basename %}
<h{{h}}>{{case}}</h{{h}}>
{{ fgsea_report(casedir, h + 1) }}
{% endfor %}
{%- endmacro -%}

{%- macro head_job(job) -%}
{% if in.casefile %}
{%  set name = in.casefile | toml_load | attr: "name" %}
{% else %}
{%  set name = envs.cases | attr: "name" %}
{% endif %}
<h1>{{name | escape}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
