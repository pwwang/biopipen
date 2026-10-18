{% from "utils/misc.liq" import report_jobs -%}
<script>
    import { Iframe } from "$libs";
</script>

{%- macro report_job(job, h=1) -%}
    {%- set qc_report = str(job.out.outdir) + "/outs/qc_report.html" -%}
    {%- if qc_report | exists -%}
    <h{{h}}>QC Report (Run-level)</h{{h}}>
    <Iframe
        src="{{qc_report}}"
        width="100%"
        frameborder="0"
        style="min-height: 60vh" />
    {%- endif -%}
    {%- set per_sample_dir = str(job.out.outdir) + "/outs/per_sample_outs" -%}
    {%- for web_summary in per_sample_dir | glob: "**" | sort -%}
        {%- if web_summary | basename == "web_summary.html" -%}
        {%- set sample_id = web_summary | replace: "/web_summary.html", "" | basename -%}
        <h{{h}}>Sample: {{sample_id | escape}}</h{{h}}>
        <Iframe
            src="{{web_summary}}"
            width="100%"
            frameborder="0"
            style="min-height: 60vh" />
        {%- endif -%}
    {%- endfor -%}
{%- endmacro -%}

{%- macro head_job(job) -%}
    <h1>{{job.out.outdir | basename | escape}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
