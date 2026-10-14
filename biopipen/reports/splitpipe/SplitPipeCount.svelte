{% from "utils/misc.liq" import report_jobs -%}
<script>
    import { Iframe } from "$libs";
</script>

{%- macro report_job(job, h=1) -%}
    {%- set report_files = job.out.outdir | glob: "**/*.html" -%}
    {%- if not report_files -%}
    <p>No HTML report was found in <code>{{job.out.outdir}}</code>.</p>
    {%- endif -%}
    {%- for report_file in report_files -%}
    <h{{h}}>{{report_file | basename | replace: "_analysis_summary.html", "" | escape}}</h{{h}}>
    <Iframe
        src="{{report_file}}"
        width="100%"
        frameborder="0"
        style="min-height: 60vh" />
    <p> </p>
    {%- endfor -%}
{%- endmacro -%}

{%- macro head_job(job) -%}
    <h1>{{job.out.outdir | basename | escape}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
