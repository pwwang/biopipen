{% from "utils/misc.liq" import report_jobs, table_of_images -%}
<script>
    import { DataTable, Image } from "@@";
</script>

{%- macro report_job(job, h=1) -%}
<h{{h}}>Summary table</h{{h}}>

<DataTable src={{job.out.outdir | joinpaths: "summary.txt" | quote}}
    data={ {{job.out.outdir | joinpaths: "summary.txt" | datatable: sep="\t" }} } />

<h{{h}}>Bar plots</h{{h}}>
{{ table_of_images(glob(job.out.outdir, "*.png")) }}

{%- endmacro -%}

{%- macro head_job(job) -%}
<h1>{{job.in.srtobj | stem}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}