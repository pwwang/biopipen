{% from "utils/misc.liq" import report_jobs, table_of_images -%}
<script>
    import { Image, DataTable } from "@@";
</script>

{%- macro report_job(job, h=1) -%}
<h{{h}}>CNV called by Control-FREEC</h{{h}}>

<p>Window = {{envs.args.general.window}}</p>

{% assign prefix = job.in.bamfile | basename %}
{% assign outprefix = job.out.outdir | joinpaths: "FREEC-output", prefix %}

<DataTable src={{ outprefix | append: "_CNVs" | quote }}
    data={ {{ outprefix | append: "_CNVs" | datatable: sep="\t", nrows=100, header=None, index_col=None }} } />

{{table_of_images(
    [outprefix + "_ratio.txt.png", outprefix + "_ratio.txt.log2.png"],
    ["Copy Number", "Copy Number (log2)"]
)}}

{%- endmacro -%}

{%- macro head_job(job) -%}
<h1>{{job.in.bamfile | stem }}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}

