{% from "utils/misc.liq" import report_jobs -%}
<script>
    import { Image } from "$libs";
    import { Tile } from "$ccs";
</script>

{%- macro report_job(job, h=1) -%}
<h{{h}}>Consistency report</h{{h}}>

<Tile>
<pre>
{{job.out.outdir | joinpaths: "consistency.txt" | read}}
</pre>
</Tile>

{%- if job.out.outdir | joinpaths: "consistency.png" | exists -%}
<h{{h}}>Heatmap of base CNV presence for each sample</h{{h}}>

<Image src={{job.out.outdir | joinpaths: "consistency.png" | quote}} />
{%- endif -%}
{%- endmacro -%}


{%- macro head_job(job) -%}
<h1>Job {{job.index}}</h1>
{%- endmacro -%}

{{ report_jobs(jobs, head_job, report_job) }}
