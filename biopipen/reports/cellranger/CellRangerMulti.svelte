{% from "utils/misc.liq" import report_jobs -%}
<script>
    import { Iframe } from "$libs";
</script>

{%- macro report_job(job, h=1) -%}
    {%- set qc_report = job.out.outdir + "/outs/qc_report.html" -%}
    {%- if qc_report | exists -%}
    <h{{h + 1}}>QC Report (Run-level)</h{{h + 1}}>
    <Iframe
        src="{{qc_report}}"
        width="100%"
        frameborder="0"
        style="min-height: 60vh" />
    {%- endif -%}
    {%- set per_sample_dir = job.out.outdir + "/outs/per_sample_outs" -%}
    {%- for sample_dir in per_sample_dir | glob: "*" | sort -%}
        {%- set sample_id = sample_dir | basename -%}
        {%- set web_summary = sample_dir + "/web_summary.html" -%}
        {%- if web_summary | exists -%}
        <h{{h + 1}}>Sample: {{sample_id | escape}}</h{{h + 1}}>
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
