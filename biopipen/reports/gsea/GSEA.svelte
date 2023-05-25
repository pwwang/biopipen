{% from "utils/misc.liq" import report_jobs -%}
{% from "utils/gsea.liq" import gsea_report -%}
<script>
    import { Image, DataTable } from "$libs";
    import { Tile } from "$ccs";
</script>

{%- macro report_job(job, h=1) -%}
{{ gsea_report(job.out.outdir, h, envs) }}
{%- endmacro -%}

{%- macro head_job(job) -%}
<h1>{{job.in.infile | stem}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
