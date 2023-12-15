{% from "utils/misc.liq" import report_jobs, table_of_images -%}

{%- macro report_job(job, h=1) -%}
    <iframe
        src="{{job.out.outdir}}/outs/web_summary.html"
        width="100%"
        frameborder="0"
        style="min-height: 60vh"></iframe>
{%- endmacro -%}

{%- macro head_job(job) -%}
    <h1>{{job.out.outdir | basename | escape}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
