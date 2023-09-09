{% from "utils/gsea.liq" import fgsea_report -%}
{% from "utils/misc.liq" import report_jobs -%}
<script>
    import { Image, DataTable } from "$libs";
</script>

{%- macro report_job(job, h=1) -%}
    {%- set secdirs = job.out.outdir | glob: "*" -%}
    {%- if len(secdirs) == 1 -%}
        {%- set secname = secdirs | first | basename -%}
        {%- for casedir in secdirs | first | glob: "*" -%}
            {%- if secname == "DEFAULT" -%}
                <h{{h}}>{{ casedir | basename | escape }}</h{{h}}>
            {%- else -%}
                <h{{h}}>{{secname | escape }} - {{ casedir | basename | escape }}</h{{h}}>
            {%- endif -%}
            {{ fgsea_report(casedir, h + 1) }}
        {%- endfor -%}
    {%- else -%}
        {%- for secdir in secdirs -%}
            {%- set sec = secdir | basename -%}
            <h{{h}}>{{sec | escape}}</h{{h}}>
            {%- for casedir in secdir | glob: "*" -%}
                <h{{h+1}}>{{casedir | basename | escape}}</h{{h+1}}>
                {{ fgsea_report(casedir, h + 2) }}
            {%- endfor -%}
        {%- endfor -%}
    {%- endif -%}
{%- endmacro -%}

{%- macro head_job(job) -%}
<h1>{{job.in.srtobj | stem0 | escape}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
