{% from "utils/misc.liq" import report_jobs -%}
<script>
    import { Image, DataTable } from "@@";
</script>

{%- macro report_job(job, h=1) -%}
<h{{h}}>TCR Clustering using {{envs.tool | upper}}</h{{h}}>

<DataTable src={{ job.out.clusterfile | quote }}
    data={ {{ job.out.clusterfile | datatable: sep="\t", nrows=100 }} } />

<h{{h+1}}>Heatmap of samples by shared TCR clusters</h{{h+1}}>

<DataTable src={{ job.outdir | joinpaths: "plotdata.csv" | quote }}
    data={ {{ job.outdir | joinpaths: "plotdata.csv" | datatable: sep="\t", index_col=0 }} } />

<Image src={{job.out.heatmap | quote}} />

{%- endmacro -%}

{%- macro head_job(job) -%}
<h1>{{job.in.immfile | stem | replace: ".immunarch", ""}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}

