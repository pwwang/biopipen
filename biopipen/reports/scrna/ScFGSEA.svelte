{% from "utils/gsea.liq" import fgsea_report -%}
{% from "utils/misc.liq" import report_jobs -%}
<script>
    import { Image, DataTable } from "$libs";
</script>

{%- macro report_job(job, h=1) -%}
{% for casedir in job.out.outdir | glob: "*" %}
    {%- set case = casedir | basename -%}
    <h{{h}}>{{case}}</h{{h}}>
    {%- if casedir | joinpaths: "percluster" | isfile -%}
        {%- for cldir in casedir | glob: "*" -%}
            {%- if basename(cldir) == "percluster" -%}
                {%- continue -%}
            {%- endif -%}
            <h{{h+1}}>{{cldir | basename}}</h{{h+1}}>
            {{ fgsea_report(cldir, h + 2) }}
        {%- endfor -%}
    {%- else -%}
        {{ fgsea_report(casedir, h + 1) }}
    {%- endif -%}
{% endfor %}
{%- endmacro -%}

{%- macro head_job(job) -%}
{% if job.in.casefile %}
{%  set name = job.in.casefile | config: "toml" | attr: "name" %}
{% else %}
{%  set name = envs | attr: "name" %}
{% endif %}
<h1>{{name | escape}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
