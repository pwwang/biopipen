{% from "utils/misc.liq" import report_jobs -%}
{% from_ os import path %}
<script>
    import { Image, DataTable, Descr } from "$libs";
</script>

{%- macro report_job(job, h=1) -%}
    {{ job | render_job: h=h }}
{%- endmacro -%}

{%- macro head_job(job) -%}
    <h1>{{job.in.metafile | stem}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
