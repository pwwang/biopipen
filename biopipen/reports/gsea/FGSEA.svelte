{% from "utils/gsea.liq" import fgsea_report -%}
{% from "utils/misc.liq" import report_jobs -%}
<script>
    import { Image, DataTable } from "@@";
</script>

{%- macro report_job(job, h=1) -%}
{{ fgsea_report(job.out.outdir, h) }}
{%- endmacro -%}

{%- macro head_job(job) -%}
<h1>{{job.out.outdir | stem | replace: ".fgsea", ""}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
